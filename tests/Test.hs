module Main (
	main
	) where

import Control.Lens
import Linear.NURBS
import Linear.V2
import Test.Hspec

-- | Simple NURBS of degree 3
test₁ ∷ NURBS V2 Double
test₁ = nurbs 3 [
	V2 0.0 0.0,
	V2 10.0 0.0,
	V2 10.0 10.0,
	V2 20.0 20.0,
	V2 0.0 20.0,
	V2 (-20.0) 0.0]

-- | Another NURBS of degree 3
test₂ ∷ NURBS V2 Double
test₂ = nurbs 3 [
	V2 (-20.0) 0.0,
	V2 (-20.0) (-20.0),
	V2 0.0 (-40.0),
	V2 20.0 20.0]

-- | Make test₁ periodic
test₀ ∷ NURBS V2 Double
test₀ = set periodic True test₁

main ∷ IO ()
main = hspec $ do
	describe "evaluate point" $ do
		it "should start from first control point" $
			(test₁ `eval` 0.0) ≃ (test₁ ^?! wpoints . _head . wpoint)
		it "should end at last control point" $
			(test₁ `eval` 1.0) ≃ (test₁ ^?! wpoints . _last . wpoint)
	describe "insert knot" $ do
		it "should not change nurbs curve" $
			insertKnots [(1, 0.1), (2, 0.3)] test₁ ≃ test₁
	describe "remove knots" $ do
		it "should not change nurbs curve" $
			removeKnots [(1, 0.1), (2, 0.3)] (insertKnots [(1, 0.1), (2, 0.3)] test₁) ≃ test₁
	describe "purge knots" $ do
		it "should not change nurbs curve" $
			purgeKnots (insertKnots [(1, 0.1), (2, 0.6)] test₁) ≃ test₁
	describe "split" $ do
		it "should work as cut" $
			snd (split 0.4 test₁) ≃ cut (Span 0.4 1.0) test₁
	describe "normalize" $ do
		it "should not affect curve" $
			cut (Span 0.2 0.8) test₁ ≃ normalizeKnot (cut (Span 0.2 0.8) test₁)
	describe "joint" $ do
		it "should joint cutted nurbs" $
			uncurry joint (split 0.3 test₁) ≃ Just test₁
		it "should cut jointed nurbs" $
			(cut (Span 0.0 1.0) <$> (test₁ ⊕ test₂)) ≃ Just test₁
		it "should cut jointed nurbs" $
			(cut (Span 1.0 2.0) <$> (test₁ ⊕ test₂)) ≃ Just test₂
	describe "periodic" $ do
		it "can be broken into simple nurbs" $
			breakLoop 0.0 test₀ ≃ test₀
		it "can be broken in any place" $
			uncurry (flip (⊕)) (split 0.5 (breakLoop 0.0 test₀)) ≃ Just (breakLoop 0.5 test₀)
