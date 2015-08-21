{-# LANGUAGE TemplateHaskell, TypeFamilies #-}

module Linear.NURBS.Types (
	Weight(..), weightPoint, weightValue, ofWeight, weight, wpoint,
	Span(..), spanStart, spanEnd,
	KnotData(..), knotData, knotDataAt, knotDataSpan,
	NURBS(..), nurbsPoints, nurbsKnot, nurbsKnoti, nurbsDegree
	) where

import Prelude.Unicode

import Control.Lens
import Linear.Vector hiding (basis)
import Linear.Affine
import Linear.Metric

-- | Point with weight
data Weight f a = Weight { _weightPoint ∷ f a, _weightValue ∷ a } deriving (Eq, Ord, Read, Show)

makeLenses ''Weight

-- | Make point with weight
ofWeight ∷ Additive f ⇒ f a → a → Weight f a
pt `ofWeight` w = Weight pt w

-- | Weight lens
weight ∷ (Additive f, Fractional a) ⇒ Lens' (Weight f a) a
weight = lens fromw tow where
	fromw (Weight _ w) = w
	tow (Weight pt w) w' = Weight ((w' / w) *^ pt) w'

-- | Point lens
wpoint ∷ (Additive f, Additive g, Fractional a) ⇒ Lens (Weight f a) (Weight g a) (f a) (g a)
wpoint = lens fromw tow where
	fromw (Weight pt w) = (1.0 / w) *^ pt
	tow (Weight _ w) pt' = Weight (w *^ pt') w

instance Functor f ⇒ Functor (Weight f) where
	fmap f (Weight pt w) = Weight (fmap f pt) (f w)	

instance Traversable f ⇒ Traversable (Weight f) where
	traverse f (Weight pt w) = Weight <$> traverse f pt <*> f w

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

-- | Knot span
data Span a = Span {
	_spanStart ∷ a,
	_spanEnd ∷ a }
		deriving (Eq, Ord, Read)

makeLenses ''Span

instance Functor Span where
	fmap f (Span s e) = Span (f s) (f e)

instance Foldable Span where
	foldMap f (Span s e) = f s `mappend` f e

instance Traversable Span where
	traverse f (Span s e) = Span <$> f s <*> f e

instance Show a ⇒ Show (Span a) where
	show (Span s e) = show (s, e)

-- | Knot evaluation data, used to compute basis functions
data KnotData a = KnotData {
	_knotDataAt ∷ a,
	_knotDataSpan ∷ Span a,
	_knotData ∷ [(Span a, a)] }
		deriving (Eq, Ord, Read, Show)

makeLenses ''KnotData

-- | NURBS
data NURBS f a = NURBS {
	_nurbsPoints ∷ [Weight f a],
	_nurbsKnot ∷ [a],
	_nurbsDegree ∷ Int }
		deriving (Eq, Ord, Read, Show)

makeLenses ''NURBS

instance Functor f ⇒ Functor (NURBS f) where
	fmap f (NURBS pts k d) = NURBS (map (fmap f) pts) (map f k) d

instance Foldable f ⇒ Foldable (NURBS f) where
	foldMap f (NURBS pts k _) = mconcat (map (foldMap f) pts) `mappend` mconcat (map f k)

nurbsKnoti ∷ Lens' (NURBS f a) [a]
nurbsKnoti = lens fromn ton where
	fromn (NURBS wpts k d)
		| length k ≡ succ (length wpts) = k
		| otherwise = drop d $ reverse $ drop d $ reverse k
	ton (NURBS wpts k d) k'
		| length k ≡ succ (length wpts) = NURBS wpts k' d
		| otherwise = NURBS wpts (replicate d (k' ^?! _head) ++ k' ++ replicate d (k' ^?! _last)) d
