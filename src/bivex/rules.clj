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
            (create-rule "k4" "demethylase" 1 0 0.5 1)
            (create-rule "k27" "methyltransferase" 0 1 1 1)
            (create-rule "k27" "demethylase" 1 0 0.5 1)
            (create-rule "k4" "turnover" 1 0 0.1 0.1)
            (create-rule "k27" "turnover" 1 0 0.1 0.1)])


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

(defn get-max-rule
  "get a rule with highest total probability and if more than one, select one at random"
  [rules]
  (let [max-prob (apply max (map get-rule-prob rules))
        max-prob-idx (map #(= max-prob %) (map get-rule-prob rules))
        max-rule (remove nil? (map #(when %2 %1) rules max-prob-idx))]
    (rand-nth max-rule)))

(defn select-rule
  "among all the applicable rules, select a rule with the highest prob. If more than one, select one at random"
  [prevnuc rules]
  (get-max-rule (get-rules-both-marks prevnuc rules)))

(defn update-rules-given-marks
  "if methyl marks are present, increase corresponding methyltransferase affinity"
  [rules prevnuc_new]
  rules
  )

(defn get-key
  [midx]
  (cond (= midx 1) :k4
        :else :k27))

(defn update-rules
  [chromtape rules]
  )


;;;;;;;;TODO

;; (defn update-rules-discourage-biv
;;   "based on the new nucleosome, discourage oppositng methyltransferase"
;;   [chromtape rules nextnuc_new]
;;   (let [x (map #(get (second nextnuc_new) %) [:k4 :k27])
;;         xtest (= (apply + (x)) 1)
;;         ])
;;   (cond xtest (filter #(and (= (:action %) "methyltransferase")
;;                             (= (:class %)
;;                                (name (get-key (+ (.indexOf x 1) 1))))) rules)
;;         :else rules))


;; (defn discourage-biv
;;   [rules givenm]
  
;;   )

