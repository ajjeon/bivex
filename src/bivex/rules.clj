(ns bivex.rules
  (:require [bivex.files :as files]))

(def default-rules-file (atom "resources/rules.csv"))

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
        new_drule (assoc (assoc drule :affinity 0.01) :abundance 0.01)
        ]
    (concat [new_drule] srule)
    ))

(defn rule-recruitment
  "mimics the recruitment by an existing methyl mark"
  [givenrules changem]
  (let [drule (into {} (filter #(and (= (:action %) "methyltransferase")
                                      (= (:class %) changem)) givenrules)) 
        srule (map #(into {} %) (filter #(or (not= (:action %) "methyltransferase")
                                              (not= (:class %) changem)) givenrules))
        new_drule (assoc (assoc drule :affinity 2) :abundance 2)
        ]
    (concat [new_drule] srule)
    ))

(defn update-rules-discourage-biv
  "based on the new nucleosome, discourage oppositng methyltransferase"
  [rules givenrules nextnuc_new]
  (let [x (map #(get (second nextnuc_new) %) [:k4 :k27])
        xtest (= (apply + x) 1)]
    (cond xtest (discourage-biv rules givenrules (name (get-key (+ (.indexOf x 1) 1))))
          :else rules)))

(defn update-rules-recruitment
  "based on the previous nucleosome, encourage recruitment by existing mark"
  [rules prevnuc_new]
  (let [updatedrules (cond (= (:k4 (second prevnuc_new)) 1)
                           (rule-recruitment rules "k4")
                           :else rules)]
    (cond (= (:k27 (second prevnuc_new)) 1)
          (rule-recruitment updatedrules "k27")
          :else rules)
    ))

(defn update-rules
  [nextnuc_new prevnuc_new]
  (let [rules (files/read-in-file @default-rules-file)
        urules (update-rules-recruitment rules prevnuc_new)] ;; after every iteration, same default rules get read in
    (update-rules-discourage-biv rules urules nextnuc_new)))

