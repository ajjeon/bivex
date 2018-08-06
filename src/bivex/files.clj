(ns bivex.files
  (:require [clojure.data.csv :as csv])
  (:require [clojure.java.io :as io])
)

(defn csv-data->maps [csv-data]
  (map zipmap
       (->> (first csv-data) ;; First row is the header
            (map keyword) ;; Drop if you want string keys instead
            repeat)
       (rest csv-data)))

(defn read-from-file [mfile]
  (csv-data->maps
   (csv/read-csv
    (io/reader mfile))))

(defn update-values [m f & args]
  (into {}
        (for [[k v] m]
          [k (cond (not (or (= k :class) (= k :action)) ) (apply f v args) 
                 :else v)])))

(defn read-in-file
  [mfile]
  (vec (map #(update-values % read-string)
            (read-from-file mfile))))

(defn read-in-chromatin
  [cfile]
  (map-indexed (fn [i v] [i v]) (read-in-file cfile)))

;; (read-in-file "resources/rules.csv")
;; (read-in-chromatin "resources/chromtape.csv")


