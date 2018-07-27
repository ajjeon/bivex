(ns bivex.core)
;  (:require [bivex.chromatin :as bc] ))

(defn foo
  "I don't do a whole lot."
  [x]
  (println x "Hello, World!"))

(defn -main
  "main function to run bivex"
  []
  ()
  )

(defn create-nucleosome
  [h, k4, k27]
  {:head h
   :k4 k4
   :k27 k27}
  )

(def nucleosome (create-nucleosome 0 "0" "0"))

;;(def chromtape (clojure.string/join "_" (repeat 10 nucleosome)))


;; (def chromtape
;;   (map-indexed (fn [i v] [i v]) [(create-nucleosome 1 "0" "0") (vec (repeat 9 nucleosome))]))

(def chromtape
  (map-indexed (fn [i v] [i v]) [(create-nucleosome 1 "0" "0")
                                 (create-nucleosome 0 "0" "0")
                                 (create-nucleosome 0 "0" "0")
                                 (create-nucleosome 0 "0" "0")
                                 (create-nucleosome 0 "0" "0")
                                 (create-nucleosome 0 "0" "0")
                                 (create-nucleosome 0 "0" "0")
                                 (create-nucleosome 0 "0" "0")
                                 (create-nucleosome 0 "0" "0")
                                 (create-nucleosome 0 "0" "0")]))

(def rule-m
  "writes K4 methyl mark"
  {:class "K4_methyltransferase"
   :status "RECRUIT"
   :left "0"
   :right "A"
   :affinity 1
   :abundance 1 })

(def rule-um
  "erases K4 methyl mark"
  {:class "K4_demethylase"
   :status "RECRUIT"
   :left "A"
   :right "0"
   :affinity 1
   :abundance 1 })


(defn expand-rule-sub
  "expand a rule when it includes wildcards"
  [ruleset]
    (if (java.lang.string/.contains (:left ruleset) "*")
      (map #(clojure.string/join(concat subrule %)) ["A","I","0"]))
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
  (nth (clojure.string/split chromtape #"_") cidx)
  )

;; work on each nucleosome. if applying rule, 
(second (nth chromtape 1))

	
(defn change-chrom
  "for changing nucleosome, change. not, copy over"
  [chrom idx]
  (let [chrom idx]
      (cond (= (:head chrom) idx)
            ( print "change"
             chrom)
            :else (print "no"))
      ))

(cond (= (:head chrom) 1)


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


