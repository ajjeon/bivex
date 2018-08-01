(ns bivex.core-test
  (:require [clojure.test :refer :all]
            [bivex.core :refer :all]))

(deftest a-test
  (testing "FIXME, I fail."
    (is (= 0 1))))








;;;;;;

;; (defn update-line
;;   [k4 k27 both genex length]
;;   (j/update-graph!
;;    "Valency Chart"
;;    {:data {:y [[k4] [k27] [both] [genex]] :x [[length][length][length][length]]}
;;     :traces [0, 1, 2, 3]}))

;; (defn evaluate-chrom
;;   [chrom_in]
;;   (let [new_chromtape (cell/apply-rule (:chromtape chrom_in))
;;         y (:k4mono (cell/check-valency new_chromtape))
;;         y2 (:k27mono (cell/check-valency new_chromtape))
;;         y3 (:biv (cell/check-valency new_chromtape))
;;         genex (cell/gene-on? y y2 y3)
;;         n (inc (:n chrom_in))
;;         ]
;;     (update-line y y2 y3 genex n) ; trace valency and gene expression
;;     (plot/plot-bar new_chromtape) ; snapshop of valency
;;     (Thread/sleep 300)
;;     {:n n :chromtape new_chromtape}))

;; (j/start-jutsu!)
;; ;(def chrom_in {:k4mono [0] :k27mono [0] :biv [0] :genex [0] :chromtape chromatin/chromtape})
;; (plot/plot-line [0] [0] [0] [0])
;; (def chrom_in {:n 0 :chromtape chromatin/chromtape})

;; (last (take 500 (iterate evaluate-chrom chrom_in)))
