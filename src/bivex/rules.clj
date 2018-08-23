(ns bivex.rules
  (:require [bivex.files :as files])
  (:require [bivex.chromatin :as chromatin]))

(def default-rules-file (atom "resources/rules.csv"))
(def new-rules-file (atom "resources/new-rules.csv"))

(defn find-rules-with-match
  "find rules with matching pattern"
  [rules key value]
  (remove nil? (map #(when (= (% key) value) %) rules)))

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

(defn get-rule-prob
  "calcualtes the total probability of a rule based on affinity and abundance"
  [r]
  (* (:affinity r) (:abundance r)))

;; (defn get-max-rule ;; OBSOLETE
;;   "get a rule with highest total probability and if more than one, select one at random"
;;   [rules]
;;   (let [max-prob (apply max (map get-rule-prob rules))
;;         max-prob-idx (map #(= max-prob %) (map get-rule-prob rules))
;;         max-rule (remove nil? (map #(when %2 %1) rules max-prob-idx))]
;;     (rand-nth max-rule)))


(defn select-weighted-rule
  [rules]
  (rand-nth (vec (concat (flatten (map #(repeat (* (* (:affinity %) (:abundance %)) 100) %) rules))))))

(defn select-rule
  "among all the applicable rules, select a rule with the highest prob. If more than one, select one at random"
  [nuc rules]
  (select-weighted-rule (get-rules-both-marks nuc rules)))

;; (defn update-rules-given-marks
;;   "if methyl marks are present, increase corresponding methyltransferase affinity"
;;   [rules prevnuc_new]
;;   rules
;;   )

(defn get-key
  [midx]
  (cond (= midx 1) :k4
        :else :k27))

(defn discourage-biv
  [rules givenrules givenm]
  (let [changem (cond (= givenm "k4") "k27" :else "k4")
        drule (into {} (filter #(and (= (:action %) "methyltransferase")
                                      (= (:class %) changem)) givenrules)) 
        srule (map #(into {} %) (filter #(or (not= (:action %) "methyltransferase")
                                              (not= (:class %) changem)) rules))
        new_drule (assoc (assoc drule :affinity 0.5) :abundance 0.5)
        ]
    (concat [new_drule] srule)
    ))


(defn update-rates
  "change feedback ratio"
  [drule srule type]
  (let [ffactor (cond (= type "methyltransferase") 2
                      :else 5)]
    (concat [(assoc (assoc drule
                           :affinity (* (:affinity drule) ffactor))
                    :abundance (* (:abundance drule) ffactor))] srule)))

(defn rule-recruitment
  "mimics the recruitment by an existing methyl mark or an empty mark" ;;;;TODO
  [givenrules changem type]
  (let [drule (into {} (filter #(and (= (:action %) type)
                                      (= (:class %) changem)) givenrules)) 
        srule (map #(into {} %) (filter #(or (not= (:action %) type)
                                             (not= (:class %) changem)) givenrules))
        ]
;    (println type changem)
;    (println givenrules)

    (cond (empty? drule) givenrules
          :else (update-rates drule srule type))))

(defn update-rules-discourage-biv
  "based on the new nucleosome, discourage oppositng methyltransferase"
  [rules givenrules nextnuc_new]
  (let [x (map #(get (second nextnuc_new) %) [:k4 :k27])
        xtest (= (apply + x) 1)]
    (cond xtest (discourage-biv rules givenrules (name (get-key (+ (.indexOf x 1) 1))))
          :else rules)))

(defn update-rules-recruitment-sub
  "adjust rates accordingly to form loop"
  [mark prevnuc_new rules]
  (cond (zero? ((keyword mark) (second prevnuc_new)))
        (rule-recruitment rules mark "demethylase")
        :else (rule-recruitment rules mark "methyltransferase")))

(defn update-rules-recruitment
  "based on the previous nucleosome, encourage recruitment by existing mark"
  [rules prevnuc_new]
  (let [updatedrules (update-rules-recruitment-sub "k4" prevnuc_new rules)]
    (update-rules-recruitment-sub "k27" prevnuc_new updatedrules)
    ))

(defn update-rules
  [new_chromtape orules nuc_h]
  (let [nextnuc_new (chromatin/find-nucleosome-with-head new_chromtape)
        prevnuc_new (nth new_chromtape (first nuc_h))
        urules (update-rules-recruitment orules prevnuc_new)] ;; after every iteration, same default rules get read in
    (update-rules-discourage-biv orules urules nextnuc_new)))
