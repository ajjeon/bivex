(ns bivex.chromatin)

(defn create-nucleosome
  [h, k4, k27]
  {:head h
   :k4 k4
   :k27 k27})

(def nucleosome (create-nucleosome 0 0 0))

(def chromtape
  (map-indexed (fn [i v] [i v]) [(create-nucleosome 0 1 0)
                                 (create-nucleosome 0 1 0)
                                 (create-nucleosome 0 1 0)
                                 (create-nucleosome 0 1 0)
                                 (create-nucleosome 0 1 0)
                                 (create-nucleosome 1 0 1)
                                 (create-nucleosome 0 0 1)
                                 (create-nucleosome 0 0 1)
                                 (create-nucleosome 0 0 1)
                                 (create-nucleosome 0 0 1)]))

;; (def chromtape
;;   (map-indexed (fn [i v] [i v]) [(create-nucleosome 0 0 0)
;;                                  (create-nucleosome 0 1 0)
;;                                  (create-nucleosome 0 0 0)
;;                                  (create-nucleosome 0 1 0)
;;                                  (create-nucleosome 0 0 1)
;;                                  (create-nucleosome 1 0 1)
;;                                  (create-nucleosome 0 1 1)
;;                                  (create-nucleosome 0 1 1)
;;                                  (create-nucleosome 0 0 0)
;;                                  (create-nucleosome 0 0 0)]))


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
