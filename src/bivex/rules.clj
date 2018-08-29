(ns bivex.rules
  (:require [bivex.files :as files])
  (:require [bivex.chromatin :as chromatin]))

(def default-rules-file (atom "resources/rules.csv"))
(def new-rules-file (atom "resources/new-rules.csv"))

(defn find-rules-with-match
  "find rules with matching pattern"
  [rules key value]
  (remove nil? (map #(when (= (% key) value) %) rules)))

(defn get-rule-prob
  "calcualtes the total probability of a rule based on affinity and abundance"
  [r]
  (* (:affinity r) (:abundance r)))

(defn get-rules-mark
  "Returns rules with matching histone mark"
  [mark nuc rules]
  (let [nucmark ((keyword mark) nuc)
        subrules (find-rules-with-match rules :class mark)]
    (find-rules-with-match subrules :left nucmark)))

(defn get-rules-both-marks
  "Returns rules with matching mark patterns given by the nucleosome"
  [nuc rules]
  (let [methyls (map name (drop 1 (keys nuc)))]
    (reduce into [] (map #(get-rules-mark % nuc rules) methyls)) ))

(defn wrand 
  "Rich Hickey's (2008) solution from ants.clj. given a vector of slice sizes, returns the index of a slice given a
  random spin of a roulette wheel with compartments proportional to
  slices."
  [slices]
  (let [total (reduce + slices)
        r (rand total)]
    (loop [i 0 sum 0]
      (if (< r (+ (slices i) sum))
        i
        (recur (inc i) (+ (slices i) sum))))))

(defn select-weighted-rule
  [rules]
  (let [wchoice (wrand (vec (map #(* (:affinity %) (:abundance %)) rules)))]
;    (println wchoice)
    (nth rules wchoice)))

(defn select-rule
  "among all the applicable rules, select a rule with the highest prob. If more than one, select one at random"
  [nuc rules]
  (let [lrules rules]
    (select-weighted-rule (get-rules-both-marks nuc lrules))))


(defn get-key
  [midx]
  (cond (= midx 1) :k4
        :else :k27))

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
(defn filter-rule
  [rules mark rtype]
  (into {} (filter #(and (= (:action %) rtype)
                         (= (:class %) mark)) rules)))

(defn rule-get-rest
  [rules rule]
  (vec (map #(cond (not= rule %) %) rules)))

(defn filter-rule-rest
  [rules exclude_rules]
  (loop [idx 0
         all_rules rules]
    (if (= idx (count exclude_rules))
      (vec (remove nil? all_rules) ) 
      (recur (inc idx) (rule-get-rest all_rules (nth exclude_rules idx))))))

(defn update-rule-rates
  [rule mfactor]
  (assoc (assoc rule :affinity (* (:affinity rule) mfactor))
          :abundance (* (:abundance rule) mfactor)))

(defn print-rule
  [rule]
  (clojure.string/join "\t" (map #(get rule %) (keys rule)))
  )

(defn print-rules
  [rules]
  (let [prules (clojure.string/join "\n" (map #(print-rule %) rules))]
  (println (clojure.string/join "\t" (map #(name %) (keys (first rules)))))
  (println prules)))

(defn feedback
  [rules mark mfactor]
  (let [pfactor mfactor
        nfactor (double (/ 1 mfactor))
        action_rule (filter-rule rules mark "methyltransferase")
        daction_rule (filter-rule rules mark "demethylase")
        maintain_rule (filter-rule rules mark "maintenance1")
        maintain2_rule (filter-rule rules mark "maintenance2")
        turnover_rule (filter-rule rules mark "turnover1")
        turnover2_rule (filter-rule rules mark "turnover2")
        action_rule_new [(update-rule-rates action_rule pfactor)] 
        daction_rule_new [(update-rule-rates daction_rule nfactor)]
        maintain_rule_new [(update-rule-rates maintain_rule pfactor)]
        maintain2_rule_new [(update-rule-rates maintain2_rule nfactor)] 
;        turnover_rule_new [(update-rule-rates turnover_rule pfactor)] 
;        turnover2_rule_new [(update-rule-rates turnover2_rule nfactor)] 
;        other_rules (filter-rule-rest rules [action_rule maintain_rule])
        ]

;    (print-rules (vec (concat [action_rule_new] [maintain_rule] other_rules)))
    (vec (concat action_rule_new daction_rule_new maintain_rule_new maintain2_rule_new [turnover_rule] [turnover2_rule]))))

(defn feedback-sub
  [rules mfactors]
  (let [m4rules (cond (not= (first mfactors) 1)
                     (feedback rules "k4" (first mfactors))
                     :else
                     (vec (filter #(= (:class %) "k4") rules)))
        m27rules (cond (not= (second mfactors) 1)
                     (feedback rules "k27" (second mfactors))
                     :else
                     (vec (filter #(= (:class %) "k27") rules)))]
;    (println mfactors)
;    (println (concat m4rules m27rules)) 
    (concat m4rules m27rules)))


(defn recruitment-based-sub
  [mark prevnuc_new]
  (cond (zero? ((keyword mark) (second prevnuc_new)))
        1
;        (feedback rules mark "negative")
        :else 2 ;(feedback rules mark "positive")
        ))

(defn recruitment-based
  [prevnuc_new]
  [(recruitment-based-sub "k4" prevnuc_new)
   (recruitment-based-sub "k27" prevnuc_new)])


(defn locus-dependent
  "locus-specific recruitment of rules"
  [nextnuc_new]
  (cond (= (:locus (second nextnuc_new)) 1) [4 1]
        (= (:locus (second nextnuc_new)) 2) [1 4]
        :else [1 1]))

(defn discourage-biv
  [nextnuc_new]
  (cond (or (= (:k4 (second nextnuc_new)) 1) (= (:k27 (second nextnuc_new)) 1))
        (cond (= (:k4 (second nextnuc_new)) 1) [1 0.5]
              (= (:k27 (second nextnuc_new)) 1) [0.5 1])
        :else [1 1]))

(defn update-rules
  [new_chromtape orules nuc_h]
  (let [nextnuc_new (chromatin/find-nucleosome-with-head new_chromtape)
        prevnuc_new (nth new_chromtape (first nuc_h))
        mfactors_r (recruitment-based prevnuc_new)
        mfactors_l (locus-dependent nextnuc_new)
        mfactors_b (discourage-biv nextnuc_new)
        mfactors [ (* (first mfactors_r) (first mfactors_l) (first mfactors_b))
                  (* (second mfactors_r) (second mfactors_l) (second mfactors_b))]] ;; after every iteration, same default rules get read in
    (feedback-sub orules mfactors)))

