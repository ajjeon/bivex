(ns bivex.rules)

(defn create-rule
  [class action left right affinity abundance]
  (hash-map :class class,
            :action action,
            :left left,
            :right right,
            :affinity affinity,
            :abundance abundance))

; TODO: JSON or tab_delim to rules
(def rules [(create-rule "k4" "methyltransferase" 0 1 1 1)
;            (create-rule "k4" "demethylase" 1 0 1 1)
            (create-rule "k27" "methyltransferase" 0 1 1 1)
 ;           (create-rule "k27" "demethylase" 1 0 1 1)
            (create-rule "k4" "turnover" 1 0 0.5 0.5)
            (create-rule "k27" "turnover" 1 0 0.5 0.5)
            (create-rule "k4" "maintenance" 1 1 1 1)
            (create-rule "k27" "maintenance" 1 1 1 1)
            (create-rule "k4" "maintenance" 0 0 1 1)
            (create-rule "k27" "maintenance" 0 0 1 1)
])


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
  [givenrules givenm]
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
  [givenrules nextnuc_new]
  (let [x (map #(get (second nextnuc_new) %) [:k4 :k27])
        xtest (= (apply + x) 1)]
    (cond xtest (discourage-biv givenrules (name (get-key (+ (.indexOf x 1) 1))))
          :else rules)))

(defn update-rules-recruitment
  "based on the new nucleosome, encourage recruitment by existing mark"
  [givenrules nextnuc_new]
  (let [updatedrules (cond (= (:k4 (second nextnuc_new)) 1) (rule-recruitment givenrules "k4")
                          :else givenrules)]
    (cond (= (:k27 (second nextnuc_new)) 1) (rule-recruitment updatedrules "k27")
                          :else updatedrules)
    ))

(defn update-rules
  [rules nextnuc_new]
  (let [urules (update-rules-discourage-biv rules nextnuc_new)]
    (update-rules-recruitment urules nextnuc_new)
    ))

