(ns bivex.cell
  (:require [bivex.chromatin :as chromatin])
  (:require [bivex.rules :as rules]))

(defn change-chrom
  "apply a rule and update the nucleosome, head leaves"
  [rule nuc]
  (assoc nuc
         (keyword (:class rule)) (:right rule)
         :head 0))

(defn drop-nuc
  [no_idx chromtape]
  (remove #(= no_idx %) chromtape))

(defn get-the-rest-idx
  "get the index of the rest of the unchanging nucleosome"
  [prevnuc_idx nextnuc_idx chromtape]
  (drop-nuc nextnuc_idx (drop-nuc prevnuc_idx (range (count chromtape)))))

(defn move-head
  "assign new head position"
  [chromtape nextnuc_idx]
  (assoc (second (nth chromtape nextnuc_idx)) :head 1))

(defn apply-rule
  [chromtape]
  (let [prevnuc_idx (chromatin/find-idx-with-head chromtape)
        nextnuc_idx (chromatin/nucleo-idx-next-head chromtape)
        prevnuc (chromatin/find-nucleosome-with-head chromtape)
        rule (rules/select-rule prevnuc)
        prevnuc_new [prevnuc_idx (change-chrom rule prevnuc)] 
        nextnuc_new [nextnuc_idx (move-head chromtape nextnuc_idx)]
        new_chromtape (concat
                       (map #(nth chromtape %)
                            (get-the-rest-idx prevnuc_idx nextnuc_idx chromtape))
                       (vector nextnuc_new)
                       (vector prevnuc_new))]
    (sort new_chromtape)))

(defn check-valency
  [new_chromtape]
  (let [k4total (vec (map #(:k4 (second %)) new_chromtape))
        k27total (vec (map #(:k27 (second %)) new_chromtape))
        biv (map #(bit-and %1 %2) k4total k27total)
        xor (map #(bit-xor %1 %2) k4total k27total)
        k4mono (map #(bit-and %1 %2) k4total xor)
        k27mono (map #(bit-and %1 %2) k27total xor)]
  {:k4mono (vector (apply + k4mono)) :k27mono (vector (apply + k27mono)) :biv (vector (apply + biv))}))

(defn gene-on?
     "calls gene expression outcome"
     [k4mono k27mono biv]
     (cond (> k4mono k27mono;(max k27mono biv)
              ) [0.5]
           :else [0.1])
     )
