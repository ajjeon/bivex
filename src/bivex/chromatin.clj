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

(defn nucleo-idx-next-head
  "Decide where to move head"
  [chromtape]
  (let [i (find-idx-with-head chromtape)]
    (cond (= i 0) (+ i 1)
          (= i (- (count chromtape) 1)) (- i 1)
          :else (+ (rand-nth [-1 1]) i) )))
