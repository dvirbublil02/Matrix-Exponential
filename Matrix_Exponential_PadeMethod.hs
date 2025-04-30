module Main where

import System.IO (hGetContents, openFile, IOMode(ReadMode))
import Numeric.LinearAlgebra hiding (forwardSubst, backSubst, luDecomp)
import Prelude hiding ((<>))
import System.Clock (Clock(Monotonic), getTime, diffTimeSpec, toNanoSecs)
import Data.List.Split (splitOn)
import qualified Data.Vector.Storable as VS
import Control.Parallel (par, pseq)
import Control.Parallel.Strategies (parMap, using, parList, rseq)
import Data.List (foldl')


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
        cols = parMap rseq (solveLU l u) (toColumns b)  -- Parallelize column solving
    in fromColumns cols

luDecomp :: Matrix Double -> (Matrix Double, Matrix Double)
luDecomp a = loop 0 (a, ident (rows a))
  where
    n = rows a
    loop k (u, l)
      | k >= n = (l, u)
      | otherwise =
          let pivot = u `atIndex` (k, k)-- list[0] = list !! 0
              -- Parallelize the update of u'
              u' = accum u (+) [((i, j), -factor * (u `atIndex` (k, j)))
                                | i <- [k+1 .. n-1],
                                  j <- [k .. n-1],
                                  let factor = u `atIndex` (i, k) / pivot] `par` u  -- Ensuring that `u` is not delayed by the parallelized part
              -- Update l' sequentially to avoid synchronization issues
              l' = accum l (+) [((i, k), u `atIndex` (i, k) / pivot) | i <- [k+1 .. n-1]]
          in u' `pseq` loop (k+1) (u', l')

solveLU :: Matrix Double -> Matrix Double -> Vector Double -> Vector Double
solveLU l u b =
    let y = forwardSubst l b
        x = backSubst u y
    in x

-- forward substitution
forwardSubst :: Matrix Double -> VS.Vector Double -> VS.Vector Double
forwardSubst l b =
    let n = VS.length b
        solve acc i =
            let sumTerms = sum [l `atIndex` (i, j) * acc VS.! j | j <- [0..(i - 1)]] -- Calculate sum of known terms
                newVal = (b VS.! i - sumTerms) / l `atIndex` (i, i) -- Calculate new value for x[i]
            in acc VS.// [(i, newVal)] -- Update accumulator with new value for x[i]

    in foldl' solve (VS.replicate n 0) [0..(n - 1)] -- Use foldl' to sequentially solve for each x[i]


-- backward substitution
backSubst :: Matrix Double -> VS.Vector Double -> VS.Vector Double
backSubst u y =
    let n = VS.length y
        solve acc i =
            let i' = n - i - 1 -- Reverse index for backward substitution
                sumTerms = sum [u `atIndex` (i', j) * acc VS.! j | j <- [(i' + 1)..(n - 1)]] -- Calculate sum of known terms
                newVal = (y VS.! i' - sumTerms) / u `atIndex` (i', i') -- Calculate new value for x[i']
            in acc VS.// [(i', newVal)] -- Update accumulator with new value for x[i']

    in foldl' solve (VS.replicate n 0) [0..(n - 1)] -- Use foldl' to sequentially solve for each x[i']


-- Padé Polynomial Computation with Early Stopping
padePolynomial :: Matrix Double -> Int -> Double -> (Matrix Double, Matrix Double)
padePolynomial a q epsilon=
    let n = rows a
        identity = ident n
        
        -- Helper function to compute terms iteratively with early stopping
        computeTerms acc1 acc2 j
            | j > q = (acc1, acc2)  -- Reached maximum iterations
            | otherwise = 
                let coefN = (fromIntegral (factorial (2*q - j)) * fromIntegral (factorial q)) /
                           (fromIntegral (factorial (2*q)) * fromIntegral (factorial j) * 
                            fromIntegral (factorial (q - j)))
                    coefD = coefN * fromIntegral ((-1)^j)

                    -- computation of terms
                    termN = scale coefN (matrixPower a j)
                    termD = scale coefD (matrixPower a j)
                    
                    -- Parallel accumulation of the terms
                    newAccN = acc1 + termN 
                    newAccD = acc2 + termD
                    -- Compute Lagrange remainder
                    remainder = computeLagrangeRemainder newAccN q
                    
                in if remainder < epsilon
                      then (newAccN, newAccD)
                      else computeTerms newAccN newAccD (j + 1)
    
    in computeTerms identity identity 1

-- Matrix Exponential Calculation
padeExp :: Matrix Double -> Double -> Matrix Double
padeExp a0 epsilon =
    let initialNorm = maximum $ map sum . toLists . tr $ cmap abs a0
        (aScaled, scaleFactor) = normalizeMatrix a0 initialNorm
        q = 8  -- The predefined value of q
        (nQ, dQ) = padePolynomial aScaled q epsilon
        expAm = solveSystem dQ nQ
    in squareMatrixNTimes expAm scaleFactor

-- Manual matrix squaring function
squareMatrixNTimes :: Matrix Double -> Int -> Matrix Double
squareMatrixNTimes mat 0 = mat
squareMatrixNTimes mat n = squareMatrixNTimes (mat <> mat) (n - 1)

-- Matrix Power Function
matrixPower :: Matrix Double -> Int -> Matrix Double
matrixPower mat 0 = ident (rows mat)  -- Identity matrix when exponent is 0
matrixPower mat n
    | even n    = let half = matrixPower mat (n `div` 2)
                  in half <> half
    | otherwise = mat <> matrixPower mat (n - 1)

-- factorial function - ex product [3, 5, 2]  -- 3 * 5 * 2 = 30
factorial :: Int -> Integer
factorial n = product [1..toInteger n]

-- Lagrange Remainder
computeLagrangeRemainder :: Matrix Double -> Int -> Double
computeLagrangeRemainder a q =
  let aNorm = maximum $ map sum . toLists . tr $ cmap abs a
      maxTerm = fromIntegral (factorial (2*q)) / 
                (product [fromIntegral (factorial j) | j <- [1..q]])
      remainder = maxTerm * (aNorm ^^ (2*q + 1) / fromIntegral (factorial (2*q + 1)))
  in remainder

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
main :: IO ()
main = do
    putStrLn "Reading matrix from TXT file..."
    a <- readTXT "exp_data.txt"
    
    let originalMatrix = a
    let epsilon = 1e-6  -- Convergence threshold

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