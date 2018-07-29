(ns bivex.rules
  (:require [cheshire.core :as cheshire])
  (:require [com.rpl.specter :as specter]))

(defrecord Rule [class action left right affinity abundance])

;(cheshire/generate-string {:foo "bar" :baz 5})

(defn create-rule
  [class action left right affinity abundance]
  (map->Rule {:class class
              :action action
              :left left
              :right right
              :affinity affinity
              :abundance abundance
              }))

(def rules [(create-rule "k4" "methyltransferase" 0 1 1 1)
            (create-rule "k4" "demethylase" 1 0 1 1)
            (create-rule "k27" "methyltransferase" 0 1 1 1)
            (create-rule "k27" "demethylase" 1 0 1 1)
            (create-rule "k4" "turnover" 1 0 0.5 0.5)
            (create-rule "k27" "turnover" 1 0 0.5 0.5)])

(def ruleset (reduce (fn [aggr {:keys [class action] :as row}]
                       (update-in aggr
                                  [class action]
                                  (fnil conj {})
                                  (dissoc row :class :action)))
                     {} rules))

(def ruleset2
  {:k4
   {:methyltransferase
    {:left 0 :right 1 :affinity 1 :abundance 1}
    :demethylase
    {:left 1 :right 0 :affinity 1 :abundance 1}
    :turnover
    {:left 1 :right 0 :affinity 0.5 :abundance 0.5}}
   :k27
   {:methyltransferase
    {:left 0 :right 1 :affinity 1 :abundance 1}
    :demethylase
    {:left 1 :right 0 :affinity 1 :abundance 1}
    :turnover
    {:left 1 :right 0 :affinity 0.5 :abundance 0.5}}})

(def map-key-walker (specter/recursive-path [akey] p [ALL (if-path [FIRST #(= % akey)] LAST [LAST p])]))

(select (map-key-walker :abundance) ruleset2)
(select (map-key-walker :aaa) {:a {:aaa 3 :b {:c {:aaa 2} :aaa 1}}})

(def TREE-VALUES
	(recursive-path [] p
	  (if-path vector?
		[ALL p]
		STAY)))
(select [TREE-VALUES] ruleset)

(defn find-all-nested
  [m k]
  (->> (tree-seq map? vals m)
       (filter map?)
       (keep k)))

(defn keys-in
  "Returns a sequence of all key paths in a given map using DFS walk."
  [m]
  (letfn [(children [node]
            (let [v (get-in m node)]
              (if (map? v)
                (map (fn [x] (conj node x)) (keys v))
                [])))
          (branch? [node] (-> (children node) seq boolean))]
    (->> (keys m)
         (map vector)
         (mapcat #(tree-seq branch? children %)))))

(defn find-highest
  [m k]
  (->> (tree-seq map? vals m)
       (filter map?)
       (keep k)))

; find rules with matching left
(defn find-rule-with-key
  "find path of keys with matching one"
  [ruleset keyid]
  (remove nil?
          (map #(when (> (.indexOf % keyid) 0) %)
               (keys-in ruleset))))

; find affnity value keypath
(def subrules
  (remove nil?
          (map #(when (> (.indexOf % :affinity) 0) %)
               (keys-in ruleset))))

(let [subrules (find-rule-with-key ruleset :affinity)])
  
(def find-matching-rule
  [ruleset matchkey matchvalue]
  ())


; return the keypath to max affin
(def max-affin
  (apply max (map #(get-in ruleset %) subrules)))

(def max-affin-idx
  (map #(= max-affin %) (map #(get-in ruleset %) subrules)))

(def max-affin-rule
  (map #(when max-affin-idx %) subrules))

; select a rule at random
(rand-nth max-affin-rule)

