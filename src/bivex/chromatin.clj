(ns bivex.chromatin
  (:require [bivex.files :as files]))

(def chromatin-file (atom "resources/chromtape.csv"))

(defn find-nucleosome-with-head*
  "Find a nucleosome with head. Returns both idx and item"
  [chromtape]
  (remove nil? (map #(when (= (:head (second %)) 1) %) chromtape)))

(defn find-nucleosome-with-head
  "Find a nucleosome with head. Returns the item"
  [chromtape]
  (second (first (find-nucleosome-with-head* chromtape)) ))

(defn find-idx-with-head
  "Find the position of the nucleosome with head"
  [chromtape]
  (first (first (find-nucleosome-with-head* chromtape))))

(defn find-idx-with-inter
  "Find the position of the nucleosome with interaction mark"
  [chromtape i]
  (let [allinter (remove nil? (map #(when (= (:inter (second %)) 1) %) chromtape))]
    (map first allinter)
    ))

(defn nucleo-idx-next-head-linear
  "Decide where to move head - linear"
  [chromtape i]
  (cond (= i 0) (+ i 1)
        (= i (- (count chromtape) 1)) (- i 1)
        :else (+ (rand-nth [-1 1]) i)))

(defn subby-idx
  [itemlist boolist]
  (remove nil? (map #(cond (not (nth boolist %)) (nth itemlist %)) (range 2))))

(defn nucleo-idx-next-head-loop
  "Decide where to move head - loop"
  [chromtape i]
  (let [interpair (remove nil? (map #(when (= (:inter (second %)) 1) %) chromtape))
        whichpair (map #(= % i) (map first interpair))]
    (first (first (subby-idx interpair whichpair)))))

(defn nucleo-idx-next-head
  [chromtape]
  (let [i (find-idx-with-head chromtape)]
    (cond (zero? (-> chromtape (nth i) second :inter))
          (nucleo-idx-next-head-linear chromtape i)
          :else (cond (zero? (rand-int 2))
                      (nucleo-idx-next-head-loop chromtape i)
                      :else
                      (nucleo-idx-next-head-linear chromtape i)))))


