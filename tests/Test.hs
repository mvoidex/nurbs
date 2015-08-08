module Main (
	main
	) where

import Prelude.Unicode

import Control.Applicative
import Control.Lens
import Data.Ratio
import Linear.NURBS
import Linear.Vector
import Linear.V2
import Test.Hspec

test ∷ NURBS V2 Double
test = nurbs 3 [
	V2 0.0 0.0,
	V2 10.0 0.0,
	V2 10.0 10.0,
	V2 20.0 20.0,
	V2 0.0 20.0,
	V2 (-20.0) 0.0]

test2 ∷ NURBS V2 Double
test2 = nurbs 3 [
	V2 (-20.0) 0.0,
	V2 (-20.0) (-20.0),
	V2 0.0 (-40.0),
	V2 20.0 20.0]

ctest ∷ NURBS V2 Double
ctest = set periodic True test

main ∷ IO ()
main = hspec $ do
	describe "insert knot" $ do
		it "should not change nurbs curve" $
			insertKnots [(1, 0.1), (2, 0.3)] test ≃ test
	describe "remove knots" $ do
		it "should not change nurbs curve" $
			removeKnots [(1, 0.1), (2, 0.3)] (insertKnots [(1, 0.1), (2, 0.3)] test) ≃ test
	describe "purge knots" $ do
		it "should not change nurbs curve" $
			purgeKnots (insertKnots [(1, 0.1), (2, 0.6)] test) ≃ test
	describe "split" $ do
		it "should work as cut" $
			snd (split 0.4 test) ≃ cut (Span 0.4 1.0) test
	describe "normalize" $ do
		it "should not affect curve" $
			cut (Span 0.2 0.8) test ≃ normalizeKnot (cut (Span 0.2 0.8) test)
	describe "joint" $ do
		it "should joint cutted nurbs" $
			uncurry joint (split 0.3 test) ≃ Just test
		it "should cut jointed nurbs" $
			(cut (Span 0.0 1.0) <$> (test ⊕ test2)) ≃ Just test
		it "should cut jointed nurbs" $
			(cut (Span 1.0 2.0) <$> (test ⊕ test2)) ≃ Just test2
	describe "periodic" $ do
		it "can be broken into simple nurbs" $
			breakLoop 0.0 ctest ≃ ctest
		it "can be broken in any place" $
			uncurry (flip (⊕)) (split 0.5 (breakLoop 0.0 ctest)) ≃ Just (breakLoop 0.5 ctest)
