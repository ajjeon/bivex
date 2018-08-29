(ns bivex.chromatin
  (:require [bivex.files :as files]))

(def chromatin-file (atom "resources/chromtape.csv"))

(defn find-nucleosome-with-head
  "Find a nucleosome with head. Returns both idx and item"
  [chromtape]
  (rand-nth (remove nil? (map #(when (= (:head (second %)) 1) %) chromtape))))

(defn find-idx-with-inter
  "Find the position of the nucleosome with interaction mark"
  [chromtape]
  (let [allinter (remove nil? (map #(when (= (:inter (second %)) 1) %) chromtape))]
    allinter))

(defn nucleo-idx-next-head-linear
  "Decide where to move head - linear"
  [chromtape i]
  (cond (= i 0) (+ i 1)
        (= i (- (count chromtape) 1)) (- i 1)
        :else (+ (rand-nth [-1 1]) i)))

(defn subby-idx
  [itemlist boolist]
  (remove nil? (map #(cond (not (nth boolist %)) (nth itemlist %)) (range (count itemlist)))))

(defn nucleo-idx-next-head-loop
  "Decide where to move head - loop"
  [chromtape i]
  (let [interpair (find-idx-with-inter chromtape)
        whichpair (map #(= % i) (map first interpair))]
    (first (rand-nth (subby-idx interpair whichpair)))))

(defn nucleo-idx-next-head
  [chromtape i]
  (cond (zero? (-> chromtape (nth i) second :inter))
        (nucleo-idx-next-head-linear chromtape i)
        :else (cond (zero? (rand-int 2))
                    (nucleo-idx-next-head-loop chromtape i)
                    :else
                    (nucleo-idx-next-head-linear chromtape i))))

(defn nucleo-get-rest
  [nuc_idx chromtape]
  (let [minus (cond (zero? nuc_idx) 0
                    :else (dec nuc_idx))
        m_range (range 0 minus)
        plus (cond (= nuc_idx (count chromtape)) (count chromtape)
                   :else (inc nuc_idx))
        p_range (range (inc plus) (inc (count chromtape)))
        all_range (vec (concat m_range p_range))]
    (map #(nth chromtape %) all_range)))


