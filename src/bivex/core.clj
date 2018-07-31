(ns bivex.core
  (:require [bivex.chromatin :as chromatin])
  (:require [bivex.rules :as rules])
  (:require [bivex.plot :as plot])
  (:require [bivex.cell :as cell])
  (:require [jutsu.core :as j])
  (:gen-class))


(defn evaluate-chrom
  [chrom_in]
  (let [new_chromtape (cell/apply-rule (:chromtape chrom_in))
        y (into (:k4mono chrom_in) (:k4mono (cell/check-valency new_chromtape)))
        y2 (into (:k27mono chrom_in) (:k27mono (cell/check-valency new_chromtape)))
        y3 (into (:biv chrom_in) (:biv (cell/check-valency new_chromtape)))
        genex (into (:genex chrom_in) (cell/gene-on? (last y) (last y2) (last y3)))
        ]
;    (plot/plot-line y y2 y3 genex) ; trace valency and gene expression
;    (plot/plot-bar new_chromtape) ; snapshot of valency
;    (Thread/sleep 100)
;    (println (clojure.string/join "__" (map plot/print-nucleosome (sort new_chromtape))))
    {:k4mono y :k27mono y2 :biv y3 :genex genex :chromtape new_chromtape}))

(j/start-jutsu!)
(def chrom_in {:k4mono [5] :k27mono [5] :biv [0] :genex [0] :chromtape chromatin/chromtape})
(last (take 100 (iterate evaluate-chrom chrom_in)))



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

;;;; cells iterate in parallel. capture snapshot after each iteration
(defn cell-valency
  [chromtape]
  (let [x (cell/check-valency chromtape)
      maxval [(apply max (flatten (vals x)))]
      xkey (filter (comp #{maxval} x) (keys x))]
;    (println x maxval xkey)
    xkey
  ))

(def cell_group (repeat 10000 chrom_in))
(def itern 100)

(let [t (pmap #(last (take itern (iterate evaluate-chrom %))) cell_group)
      x (map #(select-keys % [:k4mono :k27mono :biv :genex]) t)
      t2   (map #(map last (vals %)) x )
      k4mono (vec (map #(nth % 0) (vec t2)))
      k27mono (vec (map #(nth % 1) (vec t2)))
      biv (vec (map #(nth % 2) (vec t2)))
      genex (vec (map #(nth % 3) (vec t2)))]
  (println (count t))
;  (println (apply + k4mono) (apply + k27mono) (apply + biv) (apply + genex))
  (plot/plot-cell-all k4mono k27mono biv genex)
  (plot/plot-cell-sum k4mono k27mono biv genex) ;; TODO this looks a bit weird
  (plot/plot-bar-sum t)
;(map plot/plot-bar t); snapshop of valency
;  (plot/plot-cell-genex (map #(last (:genex %)) t) (count cell_group))
   )
