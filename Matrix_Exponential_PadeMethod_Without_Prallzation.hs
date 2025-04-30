{-# LANGUAGE FlexibleContexts #-}
module Main where

import System.IO (hGetContents, openFile, IOMode(ReadMode))
import Numeric.LinearAlgebra
import Prelude hiding ((<>))
import System.Clock (Clock(Monotonic), getTime, diffTimeSpec, toNanoSecs)
import Data.List.Split (splitOn)
import Data.List (foldl')

-- Utility Functions
factorial :: Int -> Integer
factorial n = product [1..toInteger n]

matrixPower :: Matrix Double -> Int -> Matrix Double
matrixPower mat 0 = ident (rows mat)  -- Identity matrix when exponent is 0
matrixPower mat n
    | even n    = let half = matrixPower mat (n `div` 2)
                  in half <> half
    | otherwise = mat <> matrixPower mat (n - 1)

-- Lagrange Remainder
computeLagrangeRemainder :: Matrix Double -> Int -> Bool -> Double
computeLagrangeRemainder a q isNumerator =
  let aNorm = maximum $ map sum . toLists . tr $ cmap abs a
      maxTerm = fromIntegral (factorial (2*q)) / 
                (product [fromIntegral (factorial j) | j <- [1..q]])
      remainder = maxTerm * (aNorm ^^ (2*q + 1) / fromIntegral (factorial (2*q + 1)))
  in remainder

-- Padé Polynomial Computation
padePolynomial :: Matrix Double -> Int -> Bool -> Matrix Double
padePolynomial a q isNumerator =
    let n = rows a
        identity = ident n
        terms = [ scale (coef j) (matrixPower a j) | j <- [1..q] ]
    in foldl (+) identity terms
  where
    coef j =
      let sign = if isNumerator then 1 else -1 :: Int
          numerator   = fromIntegral (factorial (2*q - j)) * fromIntegral (factorial q)
          denominator = fromIntegral (factorial (2*q)) *
                        fromIntegral (factorial j) *
                        fromIntegral (factorial (q - j))
          s           = fromIntegral (sign ^ j)
      in (numerator / denominator) * s

-- Matrix Normalization
normalizeMatrix :: Matrix Double -> Double -> (Matrix Double, Int)
normalizeMatrix mat norm
    | norm > 0.5 = let (scaledMat, m) = normalizeMatrix (scale 0.5 mat) (norm / 2)
                   in (scaledMat, m + 1)
    | otherwise = (mat, 0)

-- Padé Exponential Approximation
solveSystem :: Matrix Double -> Matrix Double -> Matrix Double
solveSystem a b =
    let (l, u) = luDecomp a
        cols = [solveLU l u col | col <- toColumns b]
    in fromColumns cols

luDecomp :: Matrix Double -> (Matrix Double, Matrix Double)
luDecomp a = loop 0 (a, ident (rows a))
  where
    n = rows a
    loop k (u, l)
      | k >= n = (l, u)
      | otherwise =
          let pivot = u `atIndex` (k, k)
              u' = accum u (+) [((i, j), -factor * (u `atIndex` (k, j)))
                                | i <- [k+1 .. n-1],
                                  j <- [k .. n-1],
                                  let factor = u `atIndex` (i, k) / pivot]
              l' = accum l (+) [((i, k), u `atIndex` (i, k) / pivot) | i <- [k+1 .. n-1]]
          in loop (k+1) (u', l')

solveLU :: Matrix Double -> Matrix Double -> Vector Double -> Vector Double
solveLU l u b = backSubst u $ forwardSubst l b

forwardSubst :: Matrix Double -> Vector Double -> Vector Double
forwardSubst l b =
    let n = size b
        updateAcc acc i = 
            let sumTerms = sum [l `atIndex` (i, j) * acc ! j | j <- [0..i-1]]
                newVal = (b ! i - sumTerms) / l `atIndex` (i, i)
            in accum acc (+) [(i, newVal)]
    in foldl' updateAcc (konst 0 n) [0..n-1]

backSubst :: Matrix Double -> Vector Double -> Vector Double
backSubst u y =
    let n = size y
        updateAcc acc i = 
            let i' = n - i - 1
                sumTerms = sum [u `atIndex` (i', j) * acc ! j | j <- [i'+1..n-1]]
                newVal = (y ! i' - sumTerms) / u `atIndex` (i', i')
            in accum acc (+) [(i', newVal)]
    in foldl' updateAcc (konst 0 n) [0..n-1]

padeExp :: Matrix Double -> Double -> Matrix Double
padeExp a0 epsilon =
  let initialNorm = maximum $ map sum . toLists . tr $ cmap abs a0
      (aScaled, scaleFactor) = normalizeMatrix a0 initialNorm
      
      q = 8  -- The predefined value of q

      nRemainder = computeLagrangeRemainder aScaled q True
      dRemainder = computeLagrangeRemainder aScaled q False
      totalRemainder = nRemainder + dRemainder

  in if totalRemainder > epsilon
         then error "Lagrange remainder exceeds the tolerance. Increase q or adjust epsilon."
         else do
            let nQ = padePolynomial aScaled q True
                dQ = padePolynomial aScaled q False
                expAm = solveSystem dQ nQ
            squareMatrixNTimes expAm scaleFactor

-- Manual matrix squaring function
squareMatrixNTimes :: Matrix Double -> Int -> Matrix Double
squareMatrixNTimes mat 0 = mat
squareMatrixNTimes mat n = squareMatrixNTimes (mat <> mat) (n - 1)

-- Timing wrapper function
timeFunction :: String -> a -> IO a
timeFunction funcName action = do
    putStrLn $ "Starting " ++ funcName ++ "..."
    startTime <- getTime Monotonic
    let result = action
    result `seq` return ()
    endTime <- getTime Monotonic
    let timeNs = toNanoSecs (diffTimeSpec endTime startTime)
    let timeSec = fromIntegral timeNs / 1e9
    putStrLn $ funcName ++ " Computation Time:"
    putStrLn $ "  Time (ns): " ++ show timeNs
    putStrLn $ "  Time (s): " ++ show timeSec
    return result

-- File Reading
readTXT :: FilePath -> IO (Matrix Double)
readTXT file = do
    handle <- openFile file ReadMode
    contents <- hGetContents handle
    let rows = lines contents
        matrixData = map (map read . splitOn " ") rows
    return $ fromLists matrixData

-- Main Function
-- Main Function
main :: IO ()
main = do
    putStrLn "Reading matrix from TXT file..."
    a <- readTXT "exp_data.txt"
    
    let originalMatrix = a

    let epsilon = 1e-6  -- Convergence threshold

    -- Time individual computational steps on the original matrix
    _ <- timeFunction "Factorial" (factorial 20)
    _ <- timeFunction "Matrix Power" (matrixPower originalMatrix 3)
    _ <- timeFunction "Lagrange Remainder" (computeLagrangeRemainder originalMatrix 8 True)
    _ <- timeFunction "Pade Polynomial" (padePolynomial originalMatrix 8 True)
    _ <- timeFunction "Matrix Normalization" 
            (normalizeMatrix originalMatrix (maximum $ map sum . toLists . tr $ cmap abs originalMatrix))
    _ <- timeFunction "LU Decomposition" (luDecomp originalMatrix)

    -- Matrix Exponential Methods Comparison
    putStrLn "\nComparing Matrix Exponential Methods:"
    
    padeExpMatrix <- timeFunction "Pade Exponential Method" (padeExp originalMatrix epsilon)
    putStrLn "Resulting matrix (Pade) - Top 5x5:"
    let expMatrixPade5x5 = take 5 . toLists $ padeExpMatrix
    print (map (take 5) expMatrixPade5x5)

    numericExpMatrix <- timeFunction "Numeric.LinearAlgebra Expm" (expm originalMatrix)
    putStrLn "Resulting matrix (Numeric.LinearAlgebra Expm) - Top 5x5:"
    let expMatrixNumeric5x5 = take 5 . toLists $ numericExpMatrix
    print (map (take 5) expMatrixNumeric5x5)