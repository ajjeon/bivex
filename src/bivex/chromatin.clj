(ns bivex.chromatin)

(defn create-nucleosome
  [h, k4, k27]
  {:head h
   :k4 k4
   :k27 k27})

(def nucleosome (create-nucleosome 0 0 0))

;; (def chromtape
;;   (map-indexed (fn [i v] [i v]) [(create-nucleosome 0 0 0)
;;                                  (create-nucleosome 0 0 0)
;;                                  (create-nucleosome 0 0 0)
;;                                  (create-nucleosome 0 0 0)
;;                                  (create-nucleosome 0 0 0)
;;                                  (create-nucleosome 1 0 0)
;;                                  (create-nucleosome 0 0 0)
;;                                  (create-nucleosome 0 0 0)
;;                                  (create-nucleosome 0 0 0)
;;                                  (create-nucleosome 0 0 0)]))

(def chromtape
  (map-indexed (fn [i v] [i v]) [(create-nucleosome 0 0 0)
                                 (create-nucleosome 0 1 0)
                                 (create-nucleosome 0 0 0)
                                 (create-nucleosome 0 1 0)
                                 (create-nucleosome 0 0 1)
                                 (create-nucleosome 1 0 1)
                                 (create-nucleosome 0 1 1)
                                 (create-nucleosome 0 1 1)
                                 (create-nucleosome 0 0 0)
                                 (create-nucleosome 0 0 0)]))


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


;; (def nucleosome {:n1 {:pos 1,
;;                        :head "0",
;;                        :k4 "0",
;;                        :k27 "0"}
;;                  :n2 {:pos 2,
;;                        :head "0",
;;                        :k4 "0",
;;                        :k27 "0"}
;;                  :n3 {:pos 3,
;;                        :head "0",
;;                        :k4 "0",
;;                        :k27 "0"}
;;                  :n4 {:pos 4,
;;                        :head "0",
;;                        :k4 "0",
;;                        :k27 "0"}
;;                  :n5 {:pos 5,
;;                        :head "0",
;;                        :k4 "0",
;;                        :k27 "0"}
;;                  :n6 {:pos 6,
;;                        :head "0",
;;                        :k4 "0",
;;                        :k27 "0"}
;;                  :n7 {:pos 7,
;;                        :head "0",
;;                        :k4 "0",
;;                        :k27 "0"}
;;                  :n8 {:pos 8,
;;                        :head "0",
;;                        :k4 "0",
;;                        :k27 "0"}
;;                  :n9 {:pos 9,
;;                        :head "0",
;;                        :k4 "0",
;;                        :k27 "0"}
;;                  :n10 {:pos 10,
;;                         :head "0",
;;                         :k4 "0",
;;                         :k27 "0"}
;;                  })
