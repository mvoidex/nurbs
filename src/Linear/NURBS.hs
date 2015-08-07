{-# LANGUAGE TypeFamilies, TemplateHaskell, RankNTypes, MultiParamTypeClasses, FlexibleInstances, DefaultSignatures, FlexibleContexts, UndecidableInstances #-}

module Linear.NURBS (
	binomial, size,
	Weight(..), ofWeight, weight, wpoint,
	Span(..), spanStart, spanEnd, spanId, grow, fall, coords, rangeSpan, mergeSpan,
	knotSpans, growSpans, fallSpans,
	KnotData(..), knotData, makeData, iterData, evalData,
	basis, rbasis,
	NURBS(..), eval, uniformKnot, degree, wpoints, points, knotVector, iknotVector, knotSpan, normalizeKnot, nurbs, wnurbs,
	insertKnot, insertKnots, appendPoint, prependPoint, split, cut, removeKnot, removeKnot_, removeKnots, purgeKnots,
	ndist, SimEq(..),
	joint, (⊕),

	pline, circle,
	) where

import Prelude.Unicode

import Control.Arrow
import Control.Lens
import Data.List
import Data.Maybe (fromMaybe)
import Linear.Vector hiding (basis)
import Linear.Affine
import Linear.Metric
import Linear.V2

binomial ∷ Integral a ⇒ a → a → a
binomial n k
	| k > n ∨ k < 0 = 0
	| otherwise = product [k + 1 .. n] `div` product [1 .. n - k]

size ∷ (Additive f, Floating a, Foldable f) ⇒ f a → a
size = sqrt ∘ sum ∘ fmap (^ (2 ∷ Integer))

data Weight f a = Weight (f a) a deriving (Eq, Ord, Read, Show)

ofWeight ∷ Additive f ⇒ f a → a → Weight f a
pt `ofWeight` w = Weight pt w

weight ∷ (Additive f, Fractional a) ⇒ Lens' (Weight f a) a
weight = lens fromw tow where
	fromw (Weight _ w) = w
	tow (Weight pt w) w' = Weight ((w' / w) *^ pt) w'

wpoint ∷ (Additive f, Additive g, Fractional a) ⇒ Lens (Weight f a) (Weight g a) (f a) (g a)
wpoint = lens fromw tow where
	fromw (Weight pt w) = (1.0 / w) *^ pt
	tow (Weight _ w) pt' = Weight (w *^ pt') w

instance Functor f ⇒ Functor (Weight f) where
	fmap f (Weight pt w) = Weight (fmap f pt) (f w)	

instance Additive f ⇒ Additive (Weight f) where
	zero = Weight zero 0
	Weight lx lw ^+^ Weight rx rw = Weight (lx ^+^ rx) (lw + rw)
	Weight lx lw ^-^ Weight rx rw = Weight (lx ^-^ rx) (lw - rw)
	lerp a (Weight lx lw) (Weight rx rw) = Weight (lerp a lx rx) (a * lw + (1 - a) * rw)
	liftU2 f (Weight lx lw) (Weight rx rw) = Weight (liftU2 f lx rx) (f lw rw)
	liftI2 f (Weight lx lw) (Weight rx rw) = Weight (liftI2 f lx rx) (f lw rw)

instance Affine f ⇒ Affine (Weight f) where
	type Diff (Weight f) = Weight (Diff f)
	Weight lx lw .-. Weight rx rw = Weight (lx .-. rx) (lw - rw)
	Weight lx lw .+^ Weight x w = Weight (lx .+^ x) (lw + w)
	Weight lx lw .-^ Weight x w = Weight (lx .-^ x) (lw - w)

instance Foldable f ⇒ Foldable (Weight f) where
	foldMap f (Weight x w) = foldMap f x `mappend` f w

instance Metric f ⇒ Metric (Weight f) where
	dot (Weight lx lw) (Weight rx rw) = dot lx rx + lw * rw	

data Span a = Span {
	_spanStart ∷ a,
	_spanEnd ∷ a }
		deriving (Eq, Ord, Read)

instance Show a ⇒ Show (Span a) where
	show (Span s e) = show (s, e)

-- | Piecewise constant function, returns 1 in span, 0 otherwise
spanId ∷ (Ord a, Num a) ⇒ a → Span a → a
spanId u (Span s e)
	| s ≤ u ∧ u ≤ e = 1
	| otherwise = 0

safeDiv ∷ (Eq a, Fractional a) ⇒ a → a → a
safeDiv _ 0.0 = 0.0
safeDiv x y = x / y

tailsOf ∷ Int → [a] → [[a]]
tailsOf n = filter ((≡ n) ∘ length) ∘ map (take n) ∘ tails

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

-- | Map value to span coordinates, span start mapped to 0, end to 1
coords ∷ (Eq a, Fractional a) ⇒ Span a → Iso' a a
coords (Span s e) = iso fromc toc where
	fromc u = (u - s) `safeDiv` (e - s)
	toc u' = (u' * (e - s)) + s

makeLenses ''Span

-- | Make span from knot vector
rangeSpan ∷ [a] → Span a
rangeSpan = uncurry Span ∘ (head &&& last)

-- | Merge spans
mergeSpan ∷ Ord a ⇒ Span a → Span a → Span a
mergeSpan l r = Span (min (_spanStart l) (_spanStart r)) (max (_spanEnd l) (_spanEnd r))

-- | Generate knot spans of degree
knotSpans ∷ Int → [a] → [Span a]
knotSpans d = map rangeSpan ∘ tailsOf (d + 2)

-- | Generate drow spans of degree
growSpans ∷ Int → [a] → [Span a]
growSpans d = knotSpans (pred d) ∘ init

-- | Generate fall spans of degree
fallSpans ∷ Int → [a] → [Span a]
fallSpans d = knotSpans (pred d) ∘ tail

-- | Knot evaluation data
data KnotData a = KnotData a [(Span a, a)] deriving (Eq, Ord, Read, Show)

knotData ∷ Lens' (KnotData a) [(Span a, a)]
knotData = lens fromk tok where
	fromk (KnotData _ d) = d
	tok (KnotData u _) = KnotData u

-- | Make initial knot data
makeData ∷ (Ord a, Num a) ⇒ [a] → a → KnotData a
makeData knot u = KnotData u $ map (id &&& spanId u) $ knotSpans 0 knot

-- | Eval basis function for next degree
iterData ∷ (Ord a, Fractional a) ⇒ KnotData a → KnotData a
iterData (KnotData u vs) = KnotData u $ zipWith mergeSpans vs (tail vs) where
	mergeSpans (ls, l) (rs, r) = (mergeSpan ls rs, grow u ls * l + fall u rs * r)

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

data NURBS f a = NURBS [Weight f a] [a] deriving (Eq, Ord, Read, Show)

instance Functor f ⇒ Functor (NURBS f) where
	fmap f (NURBS pts k) = NURBS (map (fmap f) pts) (map f k)

instance Foldable f ⇒ Foldable (NURBS f) where
	foldMap f (NURBS pts k) = mconcat (map (foldMap f) pts) `mappend` mconcat (map f k)

-- | Evaluate nurbs point
eval ∷ (Additive f, Ord a, Fractional a) ⇒ NURBS f a → a → f a
eval n t = foldr (^+^) zero [rbasis ws knot i deg t *^ pt | (i, pt) ← zip [0..] pts] where
	deg = n ^. degree
	knot = n ^. knotVector
	pts = n ^.. wpoints . each . wpoint
	ws = n ^.. wpoints . each . weight

-- | Generate knot of degree for points
uniformKnot ∷ Fractional a ⇒ Int → Int → [a]
uniformKnot deg pts = concat [
	replicate (succ deg) 0,
	[1 / fromIntegral (pts - deg) * fromIntegral i | i ← [1 .. pts - succ deg]],
	replicate (succ deg) 1]

degree ∷ Fractional a ⇒ Lens' (NURBS f a) Int
degree = lens fromn ton where
	fromn (NURBS wpts k) = length k - length wpts - 1
	ton n@(NURBS wpts _) d
		| d ≥ length wpts = n
		| d < 1 = n
		| otherwise = NURBS wpts (uniformKnot d $ length wpts)

wpoints ∷ Fractional a ⇒ Lens (NURBS f a) (NURBS g a) [Weight f a] [Weight g a]
wpoints = lens fromn ton where
	fromn (NURBS wpts _) = wpts
	ton n@(NURBS wpts k) wpts'
		| length wpts ≡ length wpts' = NURBS wpts' k
		| otherwise = NURBS wpts' (uniformKnot (view degree n) (length wpts'))

points ∷ (Additive f, Additive g, Fractional a) ⇒ Traversal (NURBS f a) (NURBS g a) (f a) (g a)
points = wpoints . each . wpoint

knotVector ∷ Eq a ⇒ Lens' (NURBS f a) [a]
knotVector = lens fromn ton where
	fromn (NURBS _ k) = k
	ton n@(NURBS wpts _) k'
		| length k' > length wpts * 2 ∨ length k' < length wpts + 2 = n
		| allSame (take (succ deg') k') ∧ allSame (take (succ deg') $ reverse k') = NURBS wpts k'
		| otherwise = n
		where
			deg' = length k' - length wpts - 1
			allSame ∷ Eq a ⇒ [a] → Bool
			allSame [] = True
			allSame (x:xs) = all (≡ x) xs

iknotVector ∷ (Eq a, Fractional a) ⇒ Lens' (NURBS f a) [a]
iknotVector = lens fromn ton where
	fromn n@(NURBS _ k) = drop (n ^. degree + 1) $ take (length k - n^. degree - 1) k
	ton n@(NURBS wpts k) k'
		| length k' + 2 * (n ^. degree + 1) > length wpts * 2 ∨ length k' + 2 * (n ^. degree + 1) < length wpts + 2 = n
		| otherwise = NURBS wpts (take (n ^. degree + 1) k ++  k' ++ take (n ^. degree + 1) (reverse k))

knotSpan ∷ (Eq a, Fractional a) ⇒ Lens' (NURBS f a) (Span a)
knotSpan = lens fromn ton where
	fromn (NURBS _ k) = rangeSpan k
	ton (NURBS wpts k) s = NURBS wpts (map (view norm') k) where
		norm' = coords (rangeSpan k) ∘ from (coords s)

normalizeKnot ∷ (Eq a, Fractional a) ⇒ NURBS f a → NURBS f a
normalizeKnot = set knotSpan (Span 0 1)

-- | Make nurbs of degree from points
nurbs ∷ (Additive f, Fractional a) ⇒ Int → [f a] → NURBS f a
nurbs deg pts = wnurbs deg (map (`ofWeight` 1) pts)

-- | Make nurbs of degree from weighted points
wnurbs ∷ (Additive f, Fractional a) ⇒ Int → [Weight f a] → NURBS f a
wnurbs deg pts = NURBS pts (uniformKnot deg (length pts))

-- | Insert knot
insertKnot ∷ (Additive f, Ord a, Fractional a) ⇒ a → NURBS f a → NURBS f a
insertKnot u n
	| u ≤ head (n ^. knotVector) ∨ u ≥ last (n ^. knotVector) = error "Invalid knot value"
	| otherwise = NURBS qs (sort $ u : (n ^. knotVector))
	where
		wpts = n ^. wpoints
		fs = map (fall u) $ fallSpans (n ^. degree) (n ^. knotVector)
		qs = [head wpts] ++ zipWith3 lerp fs wpts (tail wpts) ++ [last wpts]

-- | Insert knots
insertKnots ∷ (Additive f, Ord a, Fractional a) ⇒ [(Int, a)] → NURBS f a → NURBS f a
insertKnots iu n = foldr ($) n $ concat [replicate i (insertKnot u) | (i, u) ← iu]

-- | Append point
appendPoint ∷ (Eq a, Fractional a) ⇒ a → Weight f a → NURBS f a → NURBS f a
appendPoint knot_end pt n = NURBS
	(view wpoints n ++ [pt])
	(take (succ $ length $ view wpoints n) (view knotVector n) ++ replicate (view degree n + 1) knot_end)

-- | Prepend point
prependPoint ∷ (Eq a, Fractional a) ⇒ a → Weight f a → NURBS f a → NURBS f a
prependPoint knot_start pt n = NURBS
	(pt : view wpoints n)
	(replicate (view degree n + 1) knot_start ++ drop (view degree n) (view knotVector n))

-- | Split NURBS
split ∷ (Additive f, Ord a, Fractional a) ⇒ a → NURBS f a → (NURBS f a, NURBS f a)
split u n = (before, after) where
	n' = foldr ($) n $ replicate (view degree n - existed) (insertKnot u)
	before = NURBS (take (length bknots) $ view wpoints n') (bknots ++ replicate (view degree n' + 1) u) where
		bknots = takeWhile (< u) (view knotVector n')
	after = NURBS (drop (length (view wpoints n') - length aknots) $ view wpoints n') (replicate (view degree n' + 1) u ++ aknots) where
		aknots = dropWhile (≤ u) (view knotVector n')
	existed = length $ filter (≡ u) $ view knotVector n

-- | Cut NURBS
cut ∷ (Additive f, Ord a, Fractional a) ⇒ Span a → NURBS f a → NURBS f a
cut (Span l h) = snd ∘ split l ∘ fst ∘ split h

-- | Remove knot
removeKnot ∷ (Foldable f, Additive f, Ord a, Floating a, SimEq (NURBS f a)) ⇒ a → NURBS f a → Maybe (NURBS f a)
removeKnot u n
	| n' ≃ n = Just n'
	| otherwise = Nothing
	where
		n' = NURBS (pts ++ drop (succ $ length pts) wpts) knots'
		knots' = delete u $ n ^. knotVector
		fs = map (fall u) $ fallSpans (n ^. degree) knots'
		hs = takeWhile (> 0.0) $ 1 : [inv (1 - f) | f ← fs]
		wpts = n ^. wpoints
		pts = zipWith eval' qs hs_ where
			qs = tail (inits wpts)
			hs_ = tail (inits hs)
			eval' qs' hs' = foldr (^+^) zero $ zipWith (*^) (map h' (tails hs')) qs'
			h' [] = error "Impossible"
			h' (hi : his) = hi * product [1 - hk | hk ← his]
		inv 0.0 = 0.0
		inv x = 1.0 / x

-- | Try remove knot
removeKnot_ ∷ (Foldable f, Additive f, Ord a, Floating a, SimEq (NURBS f a)) ⇒ a → NURBS f a → NURBS f a
removeKnot_ u n = fromMaybe n $ removeKnot u n

-- | Remove knots
removeKnots ∷ (Foldable f, Additive f, Ord a, Floating a, SimEq (NURBS f a)) ⇒ [(Int, a)] → NURBS f a → NURBS f a
removeKnots iu n = foldr ($) n $ concat [replicate i (removeKnot_ u) | (i, u) ← iu]

-- | Try remove knots
purgeKnots ∷ (Foldable f, Additive f, Ord a, Floating a, SimEq (NURBS f a)) ⇒ NURBS f a → NURBS f a
purgeKnots n = foldr ($) n [removeKnot_ u | u ← n ^. iknotVector] where

ndist ∷ (Metric f, Ord a, Floating a) ⇒ f a → f a → a
ndist l r = distance l r / sqrt (max (norm l) (norm r))

class SimEq a where
	(≃) ∷ a → a → Bool
	x ≃ y = not (x ≄ y)
	(≄) ∷ a → a → Bool
	x ≄ y = not (x ≃ y)

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

joint ∷ (Ord a, Num a, Floating a, Foldable f, Metric f, SimEq (Weight f a), SimEq (NURBS f a)) ⇒ NURBS f a → NURBS f a → Maybe (NURBS f a)
joint l r
	| (l ^?! wpoints . _last) ≃ (r ^?! wpoints . _head) ∧ (l ^. degree ≡ r ^. degree) = Just $ purgeKnots $ NURBS ((l ^. wpoints) ++ (r ^. wpoints . _tail)) knots'
	| otherwise = Nothing
	where
		knots' = (l ^. knotVector . _init) ++ over mapped moveKnot (r ^.. knotVector . _tail . dropping deg' each)
		deg' = l ^. degree
		moveKnot k = k - (r ^?! knotVector . _head) + (l ^?! knotVector . _last)

(⊕) ∷ (Ord a, Num a, Floating a, Foldable f, Metric f, SimEq (Weight f a), SimEq (NURBS f a)) ⇒ NURBS f a → NURBS f a → Maybe (NURBS f a)
l ⊕ r = joint l r

pline ∷ (Additive f, Fractional a) ⇒ [f a] → NURBS f a
pline = nurbs 1

circle ∷ (Eq a, Floating a) ⇒ V2 a → a → NURBS V2 a
circle c r = over points move' $ set knotVector knot' $ wnurbs 2 [
	V2 r 0 `ofWeight` 1,
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
		knot' = [0, 0, 0, 0.25, 0.25, 0.5, 0.5, 0.75, 0.75, 1.0, 1.0, 1.0]
		move' pt = pt ^+^ c
