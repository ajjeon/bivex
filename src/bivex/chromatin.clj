(ns bivex.chromatin
  (:require [clojure.data.csv :as csv])
  (:require [clojure.java.io :as io]))

;; (defn create-nucleosome
;;   [h, k4, k27]
;;   {:head h
;;    :k4 k4
;;    :k27 k27})

;(def nucleosome (create-nucleosome 0 0 0))

(defn csv-data->maps [csv-data]
  (map zipmap
       (->> (first csv-data) ;; First row is the header
            (map keyword) ;; Drop if you want string keys instead
            repeat)
       (rest csv-data)))

(defn update-values [m f & args]
  (into {} (for [[k v] m] [k (apply f v args)])))

(def chromtape
  (map-indexed (fn [i v] [i v]) (vec (map #(update-values % read-string)
                                          (csv-data->maps
                                           (csv/read-csv
                                            (io/reader "resources/chromtape.csv"))))))) 

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
