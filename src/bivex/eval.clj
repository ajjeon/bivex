(ns bivex.eval
  (:require [bivex.chromatin :as chromatin])
  (:require [bivex.rules :as rules])
  (:require [bivex.plot :as plot])
  (:require [bivex.cell :as cell])
  (:require [bivex.files :as files])
  (:require [jutsu.core :as j])
  (:require [bivex.tools :as tools])
  (:gen-class))

(def save-chromtape (atom [[]])) ;; important for this atom to be updated and used only in single-cell runs
(def plot? (atom false))

(defn generate_chrom_in [rfile cfile trate]
  (let [chromtape (files/read-in-chromatin cfile)
        prom_chromtape (filter #(= (:position (second %)) 1) chromtape)
        valency (cell/check-valency prom_chromtape)
        k4mono (:k4mono valency)
        k27mono (:k27mono valency)
        biv (:biv valency)
        rules (files/read-in-file rfile)]
;    (println trate)
    {:k4mono k4mono
     :k27mono k27mono
     :biv biv
     :genex (cell/gene-on? (first k4mono) (first k27mono) (first biv) trate)
     :chromtape chromtape
     :trate trate
     :orules rules}))


(defn update-save-chromtape
  "saves the chromtape-state in save-chromtape atom"
  [chromtape]
  (reset! save-chromtape (merge @save-chromtape (tools/chromtape-state chromtape))))

(defn evaluate-chrom
  [chrom_in]
  (let [chromtape (tools/strip-head (:chromtape chrom_in))
        idx (first (rand-nth chromtape))
        nuc_h (nth chromtape idx)
        maxidx (first (last chromtape))
        adjnuc (rules/find-adjacent chromtape idx maxidx)
        new_rules (rules/update-rules adjnuc (nth chromtape idx) (:orules chrom_in))
        rule (rules/select-rule (second nuc_h) new_rules)
        new_chromtape (cond (= (:action rule) "turnover")
                            (rules/turnover chromtape rule nuc_h)
                            :else
                            (cell/apply-rule chromtape rule nuc_h))
        prom_chromtape (filter #(= (:position (second %)) 1) new_chromtape)
        y (into (:k4mono chrom_in) (:k4mono (cell/check-valency prom_chromtape)))
        y2 (into (:k27mono chrom_in) (:k27mono (cell/check-valency prom_chromtape)))
        y3 (into (:biv chrom_in) (:biv (cell/check-valency prom_chromtape)))
        genex (into (:genex chrom_in) (cell/gene-on? (last y) (last y2) (last y3) (:trate chrom_in)))]
    (update-save-chromtape new_chromtape)
    {:k4mono y
     :k27mono y2
     :biv y3
     :genex genex
     :chromtape (sort new_chromtape) 
     :orules (:orules chrom_in)
     :trate (:trate chrom_in)}
    ))

(defn one-plot
  [new_chrom_in]
  (plot/plot-line (:k4mono new_chrom_in)
                  (:k27mono new_chrom_in)
                  (:biv new_chrom_in)
                  (:genex new_chrom_in)) ; trace valency and gene expression
  (plot/plot-bar (:chromtape new_chrom_in))     ; snapshot of valency
  (Thread/sleep 100)
                                        ;  (println (clojure.string/join "__" (map plot/print-nucleosome (sort (:chromtape new_chrom_in)))))
  )

(defn evaluate-chrom-with-plot
  "plots updated chrom_in. only applies to single-cell iterations"
  [chrom_in]
  (let [new_chrom_in (evaluate-chrom chrom_in)]
    (when @plot? (one-plot new_chrom_in))
    new_chrom_in))

;;;; for single-cell simulations
(defn run-one
  "run one cell through iterations"
  [chrom_in itern]
   (last (take itern (iterate evaluate-chrom-with-plot chrom_in))))

(defn run-one-with-change
  [init_chrom_in beforeiter afteriter]
  (let [new-rules (files/read-in-file @rules/new-rules-file)
        beforerun (run-one init_chrom_in beforeiter)
        new_chrom_in (run-one (assoc beforerun :orules new-rules) afteriter)]
;      (reset! rules/default-rules-file @rules/new-rules-file)

    (plot/plot-nucleo-mat (drop 1 @save-chromtape)
                          (vec (remove nil?
                                       (map #(cond (= (:position (second %)) 1) (first %))
                                            (:chromtape init_chrom_in))))
                          (+ 10 (count (:chromtape init_chrom_in)))
                          (:genex new_chrom_in))))

;;;; for multiple-cell simulations

(defn run-many
  "run one cell through iterations"
  [chrom_in itern]
  (last (take itern (iterate evaluate-chrom chrom_in)))
)

(defn run-many-with-change
  [init_chrom_in beforeiter afteriter]
  (let [new-rules (files/read-in-file @rules/new-rules-file)
        beforerun (run-many init_chrom_in beforeiter)]
;      (reset! rules/default-rules-file @rules/new-rules-file)
    (run-many (assoc beforerun
               :orules new-rules) afteriter)))

(defn run-bulk
  "run multiple cells through iterations"
  [ncells beforeiter afteriter]
  (let [cell_group (map #(generate_chrom_in
                          @rules/default-rules-file @chromatin/chromatin-file %)
                        (vec (repeatedly ncells #(+ (rand-int 9) 1))))
        t (pmap #(run-many-with-change % beforeiter afteriter) cell_group)
        x (map #(select-keys % [:k4mono :k27mono :biv :genex]) t)
        t2 (map #(map last (vals %)) x)
        allgenex (map #(select-keys % [:genex]) t)
        k4mono (vec (map #(nth % 0) (vec t2)))
        k27mono (vec (map #(nth % 1) (vec t2)))
        biv (vec (map #(nth % 2) (vec t2)))
        genex (vec (map #(nth % 3) (vec t2)))]
;    (println allgenex)
;    (println (apply + k4mono) (apply + k27mono) (apply + biv) (apply + genex))
    (plot/plot-cell-all k4mono k27mono biv genex)
    (plot/plot-cell-sum k4mono k27mono biv genex ncells)
    (plot/plot-bar-sum t)
;    (println ncells allgenex beforeiter)
    (plot/when-switch-plot ncells allgenex beforeiter)
    (map #(last (:genex %)) allgenex) 
;    (plot/genex-box-plot allgenex)
    ))

