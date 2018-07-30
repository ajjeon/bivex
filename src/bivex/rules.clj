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
            (create-rule "k4" "demethylase" 1 0 1 1)
            (create-rule "k27" "methyltransferase" 0 1 1 1)
            (create-rule "k27" "demethylase" 1 0 1 1)
            (create-rule "k4" "turnover" 1 0 0.5 0.5)
            (create-rule "k27" "turnover" 1 0 0.5 0.5)])

(defn find-rules-with-match
  "find rules with matching pattern"
  [rules key value]
  (remove nil? (map #(when (= (% key) value) %) rules)))

(defn get-rules-mark
  "Returns rules with matching histone mark"
  [mark nuc]
  (let [nucmark ((keyword mark) nuc)
        subrules (find-rules-with-match rules :class mark)]
    (find-rules-with-match subrules :left nucmark)))

(defn get-rules-both-marks
  "Returns rules with matching mark patterns given by the nucleosome"
  [nuc]
  (let [methyls (map name (drop 1 (keys nuc)))]
    (reduce into [] (map #(get-rules-mark % nuc) methyls)) ))

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
  [prevnuc]
  (get-max-rule (get-rules-both-marks prevnuc)))

;; (defn expand-rule-sub
;;   "expand a rule when it includes wildcards"
;;   [ruleset]
;;     (if (java.lang.string/.contains (:left ruleset) "*")
;;       (map #(clojure.string/join(concat ruleset %)) ["A","I","0"]))
;;   )

;; (defn expand-rule
;;   "copying the head position to the expanded set of rules"
;;   [ruleset which]
;;   (let [x which]
;;     (cond
;;       (= x "left") (expand-rule-sub ruleset "1")
;;       :else (expand-rule-sub ruleset "0"))
;; ))


