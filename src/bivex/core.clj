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
;    (plot/plot-bar new_chromtape) ; snapshop of valency
;    (Thread/sleep 500)
;    (println (clojure.string/join "__" (map plot/print-nucleosome (sort new_chromtape))))
    {:k4mono y :k27mono y2 :biv y3 :genex genex :chromtape new_chromtape}))

;(j/start-jutsu!)
(def chrom_in {:k4mono [0] :k27mono [0] :biv [0] :genex [0] :chromtape chromatin/chromtape})
(last (take 500 (iterate evaluate-chrom chrom_in)))

;;;; cells iterate in parallel. capture snapshot after each iteration
(def cell_group (repeat 10 chrom_in))
(def itern 10)



(let [t (pmap #(last (take itern (iterate evaluate-chrom %))) cell_group)
      x (map #(select-keys % [:k4mono :k27mono :biv :genex]) t)
      t2   (map #(map last (vals %)) x ) ]
  (plot/plot-cell-all (vec (map #(nth % 0) (vec t2)))
           (vec (map #(nth % 1) (vec t2)))
           (vec (map #(nth % 2) (vec t2)))
           (vec (map #(nth % 3) (vec t2))))
  (plot/plot-cell-genex (map #(last (:genex %)) t) (count cell_group))

  )

;;TODO implement cell level valency call
;;; TODO implement weights

;; (defn call-valency
;;   [chromtape]
;;   (let [valency (cell/check-valency chromtape)
;;         k4mono (first (:k4mono valency))
;;         k27mono (first (:k27mono valency))
;;         biv (first (:biv valency))
;;         ]
;;     (println valency)
;;     (key (apply max-key max valency)) 
;; ))
;; ;given chromtape, call valency state
