(ns bivex.rules
  (:require [bivex.files :as files])
  (:require [bivex.chromatin :as chromatin])
  (:require [bivex.tools :as tools]))

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
    (vec (concat action_rule_new
                 daction_rule_new
                 maintain_rule_new
                 maintain2_rule_new
                 [turnover_rule] [turnover2_rule]))))

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


(defn find-inters
  [chromtape idx]
  (let [inters (vec (map first
                         (remove nil?
                                 (map
                                  #(cond (= (:inter (second %)) 1) %)
                                  chromtape))))]
    (vec (remove #(= % idx) inters))))

(defn select-dir
  [idx maxidx]
  (cond (zero? idx) "right"
        (= idx maxidx) "left"
        :else (rand-nth ["left" "right"])))

(defn get-left
  [idx maxidx]
  (remove #(< % 0) (range (- idx 1) idx))) ;; adjust neighbour proximity

(defn get-right
  [idx maxidx]
  (remove #(> % maxidx) (range (inc idx) (+ idx 2))))

(defn find-adjacent
  [chromtape idx maxidx]
  (let [dir (select-dir idx maxidx)
        adj (cond (= "left" dir) (get-left idx maxidx)
                  :else (get-right idx maxidx))
        inters (cond (= (:inter (second (nth chromtape idx))) 1)
                     (find-inters chromtape idx))
        adjidx (vec (remove nil? (concat adj inters)))]
    (map #(nth chromtape %) adjidx)))

(defn nuc-condition
  [nucs mark]
  (let [marksum (apply + (->> nucs
                              (map second)
                              (map (keyword mark))))
        middle (int (/ (count nucs) 2))]
    (cond (= marksum middle) 1
          (> marksum middle) 1.5
          (< marksum middle) 0.8)))

(defn recruitment-based
  [nucs]
  [(nuc-condition nucs "k4") (nuc-condition nucs "k27")])

(defn locus-dependent
  "locus-specific recruitment of rules"
  [nextnuc_new]
  (cond (= (:locus (second nextnuc_new)) 1) [2 1]
        (= (:locus (second nextnuc_new)) 2) [1 2]
        :else [1 1]))

(defn discourage-biv
  [nextnuc_new]
  (cond (or (= (:k4 (second nextnuc_new)) 1) (= (:k27 (second nextnuc_new)) 1))
        (cond (= (:k4 (second nextnuc_new)) 1) [1 0.1]
              (= (:k27 (second nextnuc_new)) 1) [0.1 1])
        :else [1 1]))

(defn update-rules
  [adjnucs nuc orules]
  (let [mfactors_r (recruitment-based adjnucs)
        mfactors_l (locus-dependent nuc)
        mfactors_b (discourage-biv nuc)
        mfactors [ (* (first mfactors_r) (first mfactors_l) (first mfactors_b))
                  (* (second mfactors_r) (second mfactors_l) (second mfactors_b))]]
    (feedback-sub orules mfactors)))


(defn turnover-update-chromtape 
  [chromtape mark nuc_mark_rand_idx]
  "erases mark from the chosen nucleosome"
  (let [nuc_new [nuc_mark_rand_idx
                 (assoc (second (nth chromtape nuc_mark_rand_idx)) mark 0)] 
        nuc_rest (tools/chromtape-rest chromtape [nuc_mark_rand_idx])]
    (sort (concat [nuc_new] nuc_rest))))

(defn turnover
  [chromtape rule nuc_h]
  (let [nohead_chromtape (tools/remove-head-chromtape chromtape nuc_h)
        mark (keyword (:class rule))
        nuc_mark (filter #(= (mark (second %)) 1) nohead_chromtape)
        nuc_mark_rand (cond (empty? nuc_mark) nil :else (rand-nth nuc_mark))]
    (cond (empty? nuc_mark) nohead_chromtape
            :else (turnover-update-chromtape nohead_chromtape mark (first nuc_mark_rand)))))



