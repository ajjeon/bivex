(ns bivex.eval
  (:require [bivex.chromatin :as chromatin])
  (:require [bivex.rules :as rules])
  (:require [bivex.plot :as plot])
  (:require [bivex.cell :as cell])
  (:require [bivex.files :as files])
  (:require [jutsu.core :as j])
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
     :rules rules
     :orules rules}))

(defn nucleo-state
  [nucleosome]
  "0:emtpy, 1:k4mono, 2:k27mono, 3:biv, 4:head"
  (let [k4mono (:k4 nucleosome)
        k27mono (cond (= (:k27 nucleosome) 1) 2 :else 0)
        biv (+ k4mono k27mono)
;        head (cond (= (:head nucleosome) 1) 4 :else 0)
        ]
    (max k4mono k27mono biv ;head
         )))

(defn chromtape-state
  [chrom_in]
  "assigns nucleo-state in a chromtape"
  (vec (->> chrom_in
              (:chromtape)
              (map #(nucleo-state (second %))))))


(defn update-save-chromtape
  "saves the chromtape-state in save-chromtape atom"
  [chrom_in]
  (reset! save-chromtape (merge @save-chromtape (chromtape-state chrom_in))))


;;;;;;;;;;;;;;;;;;

(defn find-inters
  [chromtape idx]
  (let [inters (vec (map first
                         (remove nil?
                                 (map
                                  #(cond (= (:inter (second %)) 1) %)
                                  chromtape))))]
    (vec (remove #(= % idx) inters))))

(defn select-dir
  [idx maxidx]
  (cond (zero? idx) "right"
        (= idx maxidx) "left"
        :else (rand-nth ["left" "right"])))

(defn get-left
  [idx maxidx]
  (remove #(< % 0) (range (- idx 3) idx))) ;; adjust neighbour proximity

(defn get-right
  [idx maxidx]
  (remove #(> % maxidx) (range (inc idx) (+ idx 4))))

(defn find-adjacent
  [chromtape idx maxidx]
  (let [dir (select-dir idx maxidx)
        adj (cond (= "left" dir) (get-left idx maxidx)
                  :else (get-right idx maxidx))
        inters (cond (= (:inter (second (nth chromtape idx))) 1)
                     (find-inters chromtape idx))
        adjidx (vec (remove nil? (concat adj inters)))]
    (map #(nth chromtape %) adjidx)))

(defn test2
  [chrom_in]
  (let [idx (first (rand-nth (:chromtape chrom_in)))
        maxidx (first (last (:chromtape chrom_in)))
        adjnuc (find-adjacent (:chromtape chrom_in) idx maxidx)
        ;;; adjust rule according to adjnuc pattern
        idx1 [(cond (zero? idx) nil :else (dec idx))] 
        idx2 [(cond (= (inc idx)  (count (:chromtape chrom_in))) nil :else (inc idx))] 
        inters 
        neighbours (vec (remove nil? (concat idx1 idx2 inters)))]
    (println idx idx1 idx2 inters neighbours)
    ))

(defn evaluate-chrom
  "applies a selected rule and update the chrom_in"
  [chrom_in]
  (let [idx (first (rand-nth (:chromtape chrom_in)))
        idx1 (vec (cond (zero? idx) nil :else (dec idx))) 
        idx2 (vec (cond (= (inc idx)  (count (:chromtape chrom_in))) nil :else (inc idx))) 
        inters (cond (= (:inter (second (nth (:chromtape chrom_in) idx))) 1)
                     (find-inters (:chromtape chrom_in) idx))
        neighbours (vec (remove nil? (concat idx1 idx2 inters))) 

;;;;;
        
; select a nucleosome at random
; check neighbouring nucleosomes
        

new_chrom_in (cell/apply-rule chrom_in)
        new_chromtape (:chromtape new_chrom_in)
        y (into (:k4mono chrom_in) (:k4mono (cell/check-valency new_chromtape)))
        y2 (into (:k27mono chrom_in) (:k27mono (cell/check-valency new_chromtape)))
        y3 (into (:biv chrom_in) (:biv (cell/check-valency new_chromtape)))
        genex (into (:genex chrom_in) (cell/gene-on? (last y) (last y2) (last y3) (:trate chrom_in)))]
    {:k4mono y
     :k27mono y2
     :biv y3
     :genex genex
     :chromtape new_chromtape
     :rules (:rules new_chrom_in)
     :orules (:orules chrom_in)
     :trate (:trate chrom_in)}))

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
    (update-save-chromtape new_chrom_in)
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
        new_chrom_in (run-one (assoc
                               (assoc beforerun :rules new-rules)
                               :orules new-rules) afteriter)]
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
    (run-many (assoc
               (assoc beforerun :rules new-rules)
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
    ;; (plot/plot-cell-all k4mono k27mono biv genex)
    ;; (plot/plot-cell-sum k4mono k27mono biv genex ncells)
    ;; (plot/plot-bar-sum t)
    ;; (println ncells allgenex beforeiter)
   ;; (plot/when-switch-plot ncells allgenex beforeiter)
    (map #(last (:genex %)) allgenex) 
;    (plot/genex-box-plot allgenex)
    ))

