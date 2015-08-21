{-# LANGUAGE TypeFamilies, RankNTypes, MultiParamTypeClasses, FlexibleInstances, DefaultSignatures, FlexibleContexts, UndecidableInstances #-}

module Linear.NURBS (
	binomial, size,
	spanId, spanEmpty, spanLength, inSpan, grow, fall, coords, rangeSpan, mergeSpan,
	knotSpans, growSpans, fallSpans,
	dataSpan, makeData, iterData, evalData,
	basis, rbasis,
	eval, uniformKnot, cycleKnot, periodic, degree, wpoints, points, knotVector, iknotVector, knotSpan, normalizeKnot, nurbs, wnurbs,
	insertKnot, insertKnots, appendPoint, prependPoint, split, cut, breakLoop, removeKnot, removeKnot_, removeKnots, purgeKnot, purgeKnots,
	ndist, SimEq(..),
	joint, (⊕),

	pline, circle,

	module Linear.NURBS.Types
	) where

import Prelude.Unicode

import Control.Arrow
import Control.Lens
import Data.List
import Data.Maybe (fromMaybe)
import Linear.Vector hiding (basis)
import Linear.Metric
import Linear.V2

import Linear.NURBS.Types

-- | Binomial coefficients
binomial ∷ Integral a ⇒ a → a → a
binomial n k
	| k > n ∨ k < 0 = 0
	| otherwise = product [k + 1 .. n] `div` product [1 .. n - k]

-- | Size of vector
size ∷ (Additive f, Floating a, Foldable f) ⇒ f a → a
size = sqrt ∘ sum ∘ fmap (^ (2 ∷ Integer))

-- | Piecewise constant function, returns 1 in span, 0 otherwise
spanId ∷ (Ord a, Num a) ⇒ a → Span a → a
spanId u s
	| u `inSpan` s = 1
	| otherwise = 0

-- | Is span empty
spanEmpty ∷ Eq a ⇒ Span a → Bool
spanEmpty (Span s e) = s ≡ e

-- | Check whether value in span
inSpan ∷ Ord a ⇒ a → Span a → Bool
x `inSpan` Span s e = x ≥ s ∧ x ≤ e

-- | Span length
spanLength ∷ Num a ⇒ Span a → a
spanLength (Span s e) = e - s

safeDiv ∷ (Eq a, Fractional a) ⇒ a → a → a
safeDiv _ 0.0 = 0.0
safeDiv x y = x / y

-- | Grow within span from 0 to 1
grow ∷ (Ord a, Fractional a) ⇒ a → Span a → a
grow u (Span l h)
	| u < l = 0
	| u > h = 1
	| l ≡ h = 1
	| otherwise = (u - l) `safeDiv` (h - l)

-- | Fall function, opposite to @grow@
fall ∷ (Ord a, Fractional a) ⇒ a → Span a → a
fall u s = 1 - grow u s

-- | Grop within span from 0 to 1, periodic
cycleGrow ∷ (Ord a, Fractional a) ⇒ a → a → Span a → a
cycleGrow per u s = grow (until (≥ (s ^. spanStart)) (+ per) u) s

-- | Fall within span from 1 to 0, periodic
cycleFall ∷ (Ord a, Fractional a) ⇒ a → a → Span a → a
cycleFall per u s = fall (until (≥ (s ^. spanStart)) (+ per) u) s

-- | Map value to span coordinates, span start mapped to 0, end to 1
coords ∷ (Eq a, Fractional a) ⇒ Span a → Iso' a a
coords (Span s e) = iso fromc toc where
	fromc u = (u - s) `safeDiv` (e - s)
	toc u' = (u' * (e - s)) + s

-- | Make span from knot vector
rangeSpan ∷ [a] → Span a
rangeSpan = uncurry Span ∘ (head &&& last)

-- | Merge spans
mergeSpan ∷ Ord a ⇒ Span a → Span a → Span a
mergeSpan l r = Span (min (_spanStart l) (_spanStart r)) (max (_spanEnd l) (_spanEnd r))

-- | Generate knot spans of degree
knotSpans ∷ Num a ⇒ Int → [a] → [Span a]
knotSpans n k = foldr ($) (zipWith Span k (tail k)) (replicate n merge') where
	merge' s = zipWith mergeSpan' s (tail s ++ [over traversed (+ _spanEnd sp) $ head s])
	mergeSpan' l r = Span (_spanStart l) (_spanEnd r)
	sp = rangeSpan k

-- | Generate drow spans of degree
growSpans ∷ Num a ⇒ Int → [a] → [Span a]
growSpans d = knotSpans (pred d) ∘ init

-- | Generate fall spans of degree
fallSpans ∷ Num a ⇒ Int → [a] → [Span a]
fallSpans d = knotSpans (pred d) ∘ tail

-- | Span of knot data
dataSpan ∷ (Eq a, Fractional a) ⇒ Lens' (KnotData a) (Span a)
dataSpan = lens fromk tok where
	fromk (KnotData _ s _) = s
	tok (KnotData u s d) s' = KnotData (scale u) s' (map (fmap scale *** scale) d) where
		scale x = x ^. coords s . from (coords s')

-- | Make initial knot data
makeData ∷ (Ord a, Num a) ⇒ [a] → a → KnotData a
makeData knot u = KnotData u (rangeSpan knot) $ map (id &&& spanId u) $ knotSpans 0 knot

-- | Eval basis function for next degree
iterData ∷ (Ord a, Fractional a) ⇒ KnotData a → KnotData a
iterData (KnotData u s vs) = KnotData u s $ zipWith mergeSpans vs (tail vs ++ [over (_1 . traversed) (+ _spanEnd s) $ head vs]) where
	mergeSpans (ls, l) (rs, r) = (mergeSpan ls rs, cycleGrow (spanLength s) u ls * l + cycleFall (spanLength s) u rs * r)

-- | Eval for n degree
evalData ∷ (Ord a, Fractional a) ⇒ Int → KnotData a → KnotData a
evalData n k = foldr ($) k (replicate n iterData)

-- | Nᵢ,ₙ — n-degree basis function for i-th control point
basis ∷ (Ord a, Fractional a) ⇒ [a] → Int → Int → a → a
basis knot i n u = evalData n (makeData knot u) ^?! knotData . ix i . _2

-- | Rᵢ,ₙ — n-degree rational basis function for i-th control point
rbasis ∷ (Ord a, Fractional a) ⇒ [a] → [a] → Int → Int → a → a
rbasis ws knot i n u = (dat ^?! knotData . ix i . _2) * (ws ^?! ix i) / sum (zipWith (*) (dat ^.. knotData . each . _2) ws) where
	dat = evalData n (makeData knot u)

-- | Evaluate nurbs point
eval ∷ (Additive f, Ord a, Fractional a) ⇒ NURBS f a → a → f a
eval n t = foldr (^+^) zero [rbasis ws knot i deg t *^ pt | (i, pt) ← zip [0..] pts] where
	deg = n ^. degree
	knot = n ^. knotVector
	pts = n ^.. wpoints . each . wpoint
	ws = n ^.. wpoints . each . weight

-- [0.0 .. 1.0] divided on n parts
unitRange ∷ Fractional a ⇒ Int → [a]
unitRange n = [fromIntegral i / fromIntegral n | i ← [0 .. n]]

-- | Generate knot of degree for points
uniformKnot ∷ Fractional a ⇒ Int → Int → [a]
uniformKnot deg pts = concat [
	replicate deg 0,
	unitRange (fromIntegral (pts - deg)),
	replicate deg 1]

-- | Generate cycle knot
cycleKnot ∷ Fractional a ⇒ Int → [a]
cycleKnot = unitRange

-- | Cut first and last same knots
cutKnot ∷ Int → [a] → [a]
cutKnot deg k = drop deg (take (length k - deg) k)

-- | Add first and last same knots
extendKnot ∷ Int → [a] → [a]
extendKnot deg k = replicate deg (k ^?! _head) ++ k ++ replicate deg (k ^?! _last)

-- | NURBS degree
degree ∷ Fractional a ⇒ Lens' (NURBS f a) Int
degree = lens fromn ton where
	fromn (NURBS _ _ d) = d
	ton n@(NURBS wpts k d) d'
		| n ^. periodic = NURBS wpts k d'
		| otherwise = NURBS wpts (extendKnot d' $ cutKnot d k) d'

-- | Is NURBS periodic
periodic ∷ Fractional a ⇒ Lens' (NURBS f a) Bool
periodic = lens fromn ton where
	fromn (NURBS wpts k _) = length k ≡ succ (length wpts)
	ton n@(NURBS wpts k d) p
		| (length k ≡ succ (length wpts)) ≡ p = n
		| p = NURBS wpts (cycleKnot (length wpts)) d
		| otherwise = NURBS wpts (uniformKnot d (length wpts)) d

-- | NURBS points with weights
wpoints ∷ Fractional a ⇒ Lens (NURBS f a) (NURBS g a) [Weight f a] [Weight g a]
wpoints = lens fromn ton where
	fromn (NURBS wpts _ _) = wpts
	ton (NURBS wpts k d) wpts'
		| length wpts ≡ length wpts' = NURBS wpts' k d
		| otherwise = NURBS wpts' (uniformKnot d (length wpts')) d

-- | NURBS points without weights
points ∷ (Additive f, Additive g, Fractional a) ⇒ Traversal (NURBS f a) (NURBS g a) (f a) (g a)
points = wpoints . each . wpoint

-- | NURBS knot vector
knotVector ∷ Eq a ⇒ Lens' (NURBS f a) [a]
knotVector = lens fromn ton where
	fromn (NURBS _ k _) = k
	ton n@(NURBS wpts _ d) k'
		| length k' ≡ succ (length wpts) = NURBS wpts k' d
		| length k' > length wpts * 2 ∨ length k' < length wpts + 2 = n
		| allSame (take (succ deg') k') ∧ allSame (take (succ deg') $ reverse k') = NURBS wpts k' deg'
		| otherwise = n
		where
			deg' = length k' - length wpts - 1
			allSame ∷ Eq a ⇒ [a] → Bool
			allSame [] = True
			allSame (x:xs) = all (≡ x) xs

iknotVector ∷ (Eq a, Fractional a) ⇒ Lens' (NURBS f a) [a]
iknotVector = lens fromn ton where
	fromn n@(NURBS _ k _)
		| n ^. periodic = k
		| otherwise = cutKnot (n ^. degree) k
	ton n@(NURBS wpts _ d) k'
		| (n ^. periodic) ∧ length k' ≡ succ (length wpts) = set knotVector k' n
		| length k' ≡ succ (length wpts - d) = set knotVector (extendKnot (n ^. degree) k') n
		| otherwise = n

-- | NURBS knot span
knotSpan ∷ (Eq a, Fractional a) ⇒ Lens' (NURBS f a) (Span a)
knotSpan = lens fromn ton where
	fromn (NURBS _ k _) = rangeSpan k
	ton (NURBS wpts k d) s = NURBS wpts (map (view norm') k) d where
		norm' = coords (rangeSpan k) ∘ from (coords s)

-- | Scale NURBS params to ∈ [0, 1]
normalizeKnot ∷ (Eq a, Fractional a) ⇒ NURBS f a → NURBS f a
normalizeKnot = set knotSpan (Span 0 1)

-- | Make nurbs of degree from points
nurbs ∷ (Additive f, Fractional a) ⇒ Int → [f a] → NURBS f a
nurbs deg pts = wnurbs deg (map (`ofWeight` 1) pts)

-- | Make nurbs of degree from weighted points
wnurbs ∷ (Additive f, Fractional a) ⇒ Int → [Weight f a] → NURBS f a
wnurbs deg pts = NURBS pts (uniformKnot deg (length pts)) deg

-- | Insert knot
-- qᵢ₊₁ = fᵢ⋅pᵢ + (1-fᵢ)⋅pᵢ₊₁
-- q₀ = p₀ (f₋₁ ≡ 0) (non periodic nurbs)
-- qₙ₊₁ = pₙ (fₙ ≡ 1) (non periodic nurbs)
-- fᵢ, pᵢ, pᵢ₊₁ ↦ qᵢ₊₁
insertKnot ∷ (Additive f, Ord a, Fractional a) ⇒ a → NURBS f a → NURBS f a
insertKnot u n
	| not (u `inSpan` (n ^. knotSpan)) = error "Invalid knot value"
	| otherwise = NURBS qs (sort $ u : (n ^. knotVector)) (n ^. degree)
	where
		wpts = n ^. wpoints
		fs = map (cycleFall (spanLength (n ^. knotSpan)) u) $ take (length (n ^. wpoints)) $ knotSpans (pred (n ^. degree)) (n ^. knotVector)
		-- find place to insert additional point
		-- after last non null fall param
		(growth, restart) = over (both . mapped . _1) snd $ span (uncurry (≤) ∘ fst) $ zip (zip (0.0 : fs) fs) (zip (last wpts : wpts) wpts)
		qs = map (\(f, (prev, cur)) → lerp f prev cur) $ growth ++ (set _1 0.0 (last growth)) : restart

-- | Insert knots
insertKnots ∷ (Additive f, Ord a, Fractional a) ⇒ [(Int, a)] → NURBS f a → NURBS f a
insertKnots iu n = foldr ($) n $ concat [replicate i (insertKnot u) | (i, u) ← iu]

-- | Append point
appendPoint ∷ (Eq a, Fractional a) ⇒ a → Weight f a → NURBS f a → NURBS f a
appendPoint knot_end pt =
	over nurbsPoints (++ [pt]) ∘
	over nurbsKnoti (++ [knot_end])

-- | Prepend point
prependPoint ∷ (Eq a, Fractional a) ⇒ a → Weight f a → NURBS f a → NURBS f a
prependPoint knot_start pt =
	over nurbsPoints (pt :) ∘
	over nurbsKnoti (knot_start :)

-- | Split NURBS
split ∷ (Additive f, Ord a, Fractional a) ⇒ a → NURBS f a → (NURBS f a, NURBS f a)
split u n
	| n ^. periodic = error "Can't split periodic nurbs"
	| otherwise = (before, after)
	where
		n' = insertKnots [(n ^. degree - existed, u)] n
		before = NURBS (take (length bknots) $ view wpoints n') (bknots ++ replicate (view degree n' + 1) u) (n ^. degree) where
			bknots = takeWhile (< u) (view knotVector n')
		after = NURBS (drop (length (view wpoints n') - length aknots) $ view wpoints n') (replicate (view degree n' + 1) u ++ aknots) (n ^. degree) where
			aknots = dropWhile (≤ u) (view knotVector n')
		existed = length $ filter (≡ u) $ view knotVector n

-- | Cut NURBS
cut ∷ (Additive f, Ord a, Fractional a) ⇒ Span a → NURBS f a → NURBS f a
cut (Span l h) n
	| not (n ^. periodic) = snd ∘ split l ∘ fst ∘ split h $ n
	| otherwise = fst ∘ split h ∘ breakLoop l $ n

-- | Break periodic NURBS at param
breakLoop ∷ (Additive f, Ord a, Fractional a) ⇒ a → NURBS f a → NURBS f a
breakLoop u n
	| not (n ^. periodic) = error "Can't break not periodic nurbs"
	| otherwise = NURBS (wpts' ++ [head wpts']) knot' (n ^. degree)
	where
		n' = insertKnots [(n ^. degree - existed, u)] n
		existed = length $ filter (≡ u) $ n ^. knotVector
		(knot_tail, knot_init) = break (≡ u) $ n' ^. knotVector . _init
		knot' =
			extendKnot (n' ^. degree) $
				drop (pred $ n' ^. degree) knot_init ++
				over each (+ (n' ^. knotSpan . spanEnd)) (knot_tail ++ [u])
		wpts' = rotate (pred $ length knot_tail) (n' ^. wpoints)

-- | Remove knot
-- pᵢ₊₁ = (qᵢ₊₁ - fᵢ⋅pᵢ)/(1-fᵢ) = hᵢ⋅qᵢ₊₁ + (1-hᵢ)⋅pᵢ, where hᵢ = 1/(1-fᵢ) ∧ fᵢ ≢ 1
-- if fᵢ = 1 then pᵢ₊₁ = qᵢ₊₁
-- p₀ = q₀ (h₋₁ ≡ 1)
-- pₙ = qₙ₊₁ (hₙ ≡ ∞, fₙ ≡ 1)
-- hᵢ, qᵢ₊₁, pᵢ ↦ pᵢ₊₁
removeKnot ∷ (Foldable f, Additive f, Ord a, Floating a, SimEq (NURBS f a)) ⇒ a → NURBS f a → Maybe (NURBS f a)
removeKnot u n
	| n ^. periodic = Nothing
	| n' ≃ n = Just n'
	| otherwise = Nothing
	where
		n' = NURBS (pts ++ drop (succ $ length pts) wpts) knots' (n ^. degree)
		knots' = delete u $ n ^. knotVector
		fs = map (fall u) $ fallSpans (n ^. degree) knots'
		hs = takeWhile (> 0.0) [inv (1 - f) | f ← fs]
		wpts = n ^. wpoints
		pts = head wpts : zipWith3 lerp hs (tail wpts) pts
		inv 0.0 = 0.0
		inv x = 1.0 / x

-- | Try remove knot
removeKnot_ ∷ (Foldable f, Additive f, Ord a, Floating a, SimEq (NURBS f a)) ⇒ a → NURBS f a → NURBS f a
removeKnot_ u n = fromMaybe n $ removeKnot u n

-- | Remove knots
removeKnots ∷ (Foldable f, Additive f, Ord a, Floating a, SimEq (NURBS f a)) ⇒ [(Int, a)] → NURBS f a → NURBS f a
removeKnots iu n = foldr ($) n $ concat [replicate i (removeKnot_ u) | (i, u) ← iu]

-- | Try remove knot as much times as possible
purgeKnot ∷ (Foldable f, Additive f, Ord a, Floating a, SimEq (NURBS f a)) ⇒ a → NURBS f a → NURBS f a
purgeKnot u n = fromMaybe n (purgeKnot u <$> removeKnot u n)

-- | Try remove knots
purgeKnots ∷ (Foldable f, Additive f, Ord a, Floating a, SimEq (NURBS f a)) ⇒ NURBS f a → NURBS f a
purgeKnots n = foldr ($) n [removeKnot_ u | u ← vs] where
	vs
		| n ^. periodic = n ^. knotVector
		| otherwise = cutKnot (n ^. degree + 1) $ n ^. knotVector

-- | Distance between points
ndist ∷ (Metric f, Ord a, Floating a) ⇒ f a → f a → a
ndist l r = distance l r / sqrt (max (norm l) (norm r))

class SimEq a where
	(≃) ∷ a → a → Bool
	x ≃ y = not (x ≄ y)
	(≄) ∷ a → a → Bool
	x ≄ y = not (x ≃ y)

simEq ∷ SimEq a ⇒ a → a → Bool
simEq = (≃)

simNeq ∷ SimEq a ⇒ a → a → Bool
simNeq = (≄)

instance SimEq a ⇒ SimEq (Maybe a) where
	Just x ≃ Just y = x ≃ y
	_ ≃ _ = False

instance {-# OVERLAPPABLE #-} (Ord a, Floating a) ⇒ SimEq a where
	x ≃ y = abs (x - y) ≤ 1e-6 * max x y

instance {-# OVERLAPPABLE #-} (Metric f, Ord a, Floating a, SimEq a) ⇒ SimEq (f a) where
	x ≃ y = distance x y ≤ 1e-6 * max (norm x) (norm y)

instance {-# OVERLAPS #-} (Metric f, Ord a, Floating a, SimEq a) ⇒ SimEq (Weight f a) where
	x ≃ y = distance x y ≤ 1e-6 * max (norm x) (norm y)

instance (Metric f, Ord a, Floating a, SimEq (f a)) ⇒ SimEq (NURBS f a) where
	x ≃ y = dt ≤ 1e-6 * max (norm' x) (norm' y) where
		dt = sum $ map norm $ zipWith (^-^)
			(map (eval x) $ u ^.. each . from (coords (x ^. knotSpan)))
			(map (eval y) $ u ^.. each . from (coords (y ^. knotSpan)))
		u = map ((/ fromIntegral i) ∘ fromIntegral) [0 .. i]
		i = max (length (x ^. wpoints)) (length (y ^. wpoints)) * 4
		norm' n = maximum $ map norm (n ^. wpoints)

-- | Try to joint two NURBS
joint ∷ (Ord a, Num a, Floating a, Foldable f, Metric f, SimEq (Weight f a), SimEq (NURBS f a)) ⇒ NURBS f a → NURBS f a → Maybe (NURBS f a)
joint l r
	| (l ^?! wpoints . _last) ≃ (r ^?! wpoints . _head) ∧ (l ^. degree ≡ r ^. degree) = Just $ purgeKnot (l ^?! knotVector . _last) $ NURBS ((l ^. wpoints) ++ (r ^. wpoints . _tail)) knots' (l ^. degree)
	| otherwise = Nothing
	where
		knots' = (l ^. knotVector . _init) ++ over mapped moveKnot (r ^.. knotVector . _tail . dropping deg' each)
		deg' = l ^. degree
		moveKnot k = k - (r ^?! knotVector . _head) + (l ^?! knotVector . _last)

-- | Joint
(⊕) ∷ (Ord a, Num a, Floating a, Foldable f, Metric f, SimEq (Weight f a), SimEq (NURBS f a)) ⇒ NURBS f a → NURBS f a → Maybe (NURBS f a)
l ⊕ r = joint l r

-- | Make pline NURBS
pline ∷ (Additive f, Fractional a) ⇒ [f a] → NURBS f a
pline = nurbs 1

-- | Make circle NURBS
circle ∷ (Eq a, Floating a) ⇒ V2 a → a → NURBS V2 a
circle c r = over points move' $ set nurbsKnot knot' $ wnurbs 2 [
	V2 r r `ofWeight` sq,
	V2 0 r `ofWeight` 1,
	V2 (-r) r `ofWeight` sq,
	V2 (-r) 0 `ofWeight` 1,
	V2 (-r) (-r) `ofWeight` sq,
	V2 0 (-r) `ofWeight` 1,
	V2 r (-r) `ofWeight` sq,
	V2 r 0 `ofWeight` 1]
	where
		sq = sqrt 2.0 / 2.0
		knot' = [0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0]
		move' pt = pt ^+^ c

rotate ∷ Int → [a] → [a]
rotate n l = uncurry (flip (++)) $ splitAt n' l where
	n' = until (≥ 0) (+ length l) n
