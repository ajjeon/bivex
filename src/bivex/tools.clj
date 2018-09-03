(ns bivex.tools)

(defn drop-nuc
  [no_idx chromtape_idx]
  (remove #(= no_idx %) chromtape_idx))

(defn drop-nuc-many
  [drop_idx chromtape_idx]
  (loop [idx 0
         all_idx chromtape_idx]  
    (if (= idx (count drop_idx))
      all_idx
      (recur (inc idx) (drop-nuc (nth drop_idx idx) all_idx)))))

(defn chromtape-rest
  [chromtape drop_idx]
  (let [chromtape_idx (range (count chromtape))]
    (map #(nth chromtape %) (drop-nuc-many drop_idx chromtape_idx))))

(defn cast-head
  "assign new head position"
  [chromtape nextnuc_idx]
  (assoc (second (nth chromtape nextnuc_idx)) :head 1))

(defn remove-head
  "remove head for a given nucleosome"
  [nuc]
  (assoc nuc :head 0))

(defn remove-head-chromtape
  [chromtape nuc_h]
  (let [nuc_h_new [(first nuc_h) (remove-head (second nuc_h))]
        nuc_rest (chromtape-rest chromtape [(first nuc_h)])]
    (sort (concat [nuc_h_new] nuc_rest))))

(defn change-chrom
  "apply a rule and update the nucleosome, head leaves"
  [rule nuc]
  (assoc nuc
         (keyword (:class rule)) (:right rule)
         :head 1))

(defn strip-head-nuc
  [nuc]
  [(first nuc) (assoc (second nuc) :head 0)])

(defn strip-head
  [chromtape]
  (map strip-head-nuc chromtape))

;; (defn put-head
;;   [chromtape idx]
;;   (let [others (strip-head (chromtape-rest chromtape [idx])) 
;;         newnuc  [idx (assoc (second (nth chromtape idx)) :head 1)]]
;;     (sort (concat others [newnuc]))))

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
  [chromtape]
  "assigns nucleo-state in a chromtape"
  (vec (->> chromtape
            (map #(nucleo-state (second %))))))

