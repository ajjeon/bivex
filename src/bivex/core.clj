(ns bivex.core
  (:require [bivex.rules]))
;  (:require [bivex.chromatin :as bc] ))

(defn -main
  "main function to run bivex"
  []
  ()
  )

(defn create-nucleosome
  [h, k4, k27]
  {:head h
   :k4 k4
   :k27 k27})

(def nucleosome (create-nucleosome 0 0 0))

;;(def chromtape (clojure.string/join "_" (repeat 10 nucleosome)))


;; (def chromtape
;;   (map-indexed (fn [i v] [i v]) [(create-nucleosome 1 "0" "0") (vec (repeat 9 nucleosome))]))

(def chromtape
  (map-indexed (fn [i v] [i v]) [(create-nucleosome 1 0 0)
                                 (create-nucleosome 0 0 0)
                                 (create-nucleosome 0 0 0)
                                 (create-nucleosome 0 0 0)
                                 (create-nucleosome 0 0 0)
                                 (create-nucleosome 0 0 0)
                                 (create-nucleosome 0 0 0)
                                 (create-nucleosome 0 0 0)
                                 (create-nucleosome 0 0 0)
                                 (create-nucleosome 0 0 0)]))


(defn expand-rule-sub
  "expand a rule when it includes wildcards"
  [ruleset]
    (if (java.lang.string/.contains (:left ruleset) "*")
      (map #(clojure.string/join(concat ruleset %)) ["A","I","0"]))
  )

(defn expand-rule
  "copying the head position to the expanded set of rules"
  [ruleset which]
  (let [x which]
    (cond
      (= x "left") (expand-rule-sub ruleset "1")
      :else (expand-rule-sub ruleset "0"))
))

(defn copy-chrom
  "copy the marks from input chromatin"
  [chromatin cidx]
;  (nth (clojure.string/split chromtape #"_") cidx)
  )

;; work on each nucleosome. if applying rule, 

(defn find-rule
  "find all rules with matching left"
  []
  ())

(defn select-rule
  "select a rule among all the matching rules. This steps considers the affinity and abundance"
  []
  ())



(defn change-chrom
  "for changing nucleosome, change. not, copy over"
  [chrom]
  (let [chrom chrom]
      (cond (= (:head chrom) 1)


            (print "change")
            :else (print "no"))))

(map #(change-chrom (second (nth chromtape %)))
     (take 10 (range)))


;; (if idx == changing index, change the nucleosome
;;    not, copy the nucleosome

(doseq [[idx item] chromtape]
;  (println item)
  (cond (= (:head item) "1") (change-chrom item idx)
  )
)


;; (clojure.walk/prewalk-replace [(second (nth chromtape 2)) [:head "1"]]) ;; updating values

;; (map #(vector % (inc (% chromtape))) [:a :b])

;; (vector 1 (update-in (second (nth chromtape 1)) [:head] inc))

;;;;; at each iteration, recruit-change-move-eval


(cond (= state "RECRUIT") () ;; recruit matching modifiers
      (= state "APPLY") () ;; apply changes
      (= state "MOVE") () ;; move head
      (= state "EVALUATE") () ;; evaluate - iteration, methyl proportion, gene exp outcome
      )


