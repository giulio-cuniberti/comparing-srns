## [1] Example 4.4 (Reversible reaction) ###########################################################

## ordered species: 
##   1. S
##   2. P

## ordered reactions:
##   1. S --> P
##   2. P --> S


##           (S,P)
S1 <- rbind(c(1,0),
            c(0,1))

##           (S,P)
P1 <- rbind(c(0,1),
            c(1,0))

##            (S,P)
M1 <- rbind(c(-1,0))


check1 <- checkOrdering(source.complexes = S1, product.complexes = P1, preorder.matrix = M1)

orderings1 <- findOrderings(source.complexes = S1, product.complexes = P1)



## [2] Example 4.5 (SIS) ###########################################################################

## ordered species: 
##   1. S
##   2. I

## ordered reactions:
##   1. S+I --> 2I
##   2.   I --> S


##           (S,I)
S2 <- rbind(c(1,1),
            c(0,1))

##           (S,I)
P2 <- rbind(c(0,2),
            c(1,0))


orderings2 <- findOrderings(source.complexes = S2, product.complexes = P2)



## [3] Example 4.5 (SIR) ###########################################################################

## ordered species: 
##   1. S
##   2. I
##   3. R

## ordered reactions:
##   1. S+I --> 2I
##   2.   I --> R


##           (S,I,R)
S3 <- rbind(c(1,1,0),
            c(0,1,0))

##           (S,I,R)
P3 <- rbind(c(0,2,0),
            c(0,0,1))

##            (S,I,R)
M3 <- rbind(c(-1,0,0),
             c(0,1,0))


check3 <- checkOrdering(source.complexes = S3, product.complexes = P3, preorder.matrix = M3)

orderings3 <- findOrderings(source.complexes = S3, product.complexes = P3)



## [4] Example 4.6 (Michaelis-Menten) ##############################################################

## ordered species: 
##   1. S
##   2. E
##   3. C
##   4. P

## ordered reactions:
##   1. S+E --> C
##   2.   C --> S+E
##   3.   C --> E+P


##           (S,E,C,P)
S4 <- rbind(c(1,1,0,0),
            c(0,0,1,0),
            c(0,0,1,0))

##           (S,E,C,P)
P4 <- rbind(c(0,0,1,0),
            c(1,1,0,0),
            c(0,1,0,1))


orderings4 <- findOrderings(source.complexes = S4, product.complexes = P4)



## [5] Example 4.7 (Reversible Michaelis-Menten) ###################################################

## ordered species: 
##   1. S
##   2. E
##   3. C
##   4. P

## ordered reactions:
##   1. S+E --> C
##   2.   C --> S+E
##   3.   C --> E+P
##   4.   P --> S


##           (S,E,C,P)
S5 <- rbind(c(1,1,0,0),
            c(0,0,1,0),
            c(0,0,1,0),
            c(0,0,0,1))

##           (S,E,C,P)
P5 <- rbind(c(0,0,1,0),
            c(1,1,0,0),
            c(0,1,0,1),
            c(1,0,0,0))

##            (S,E,C,P)
M5 <- rbind(c(-1,0,0,0),
             c(0,0,0,1))


check5 <- checkOrdering(source.complexes = S5, product.complexes = P5, preorder.matrix = M5)

orderings5 <- findOrderings(source.complexes = S5, product.complexes = P5)



## [6] Example 4.8 (Signaling cascade) #############################################################

## ordered species: 
##   1. P0
##   2. S1
##   3. C1
##   4. P1
##   5. S2
##   6. C2
##   7. P2
##   8. S3
##   9. C3
##  10. P3

## ordered reactions:
##   1. P0+S1 --> C1
##   2.    C1 --> P0+S1
##   3.    C1 --> P0+P1
##   4. P1+S2 --> C2
##   5.    C2 --> P1+S2
##   6.    C2 --> P1+P2
##   7. P2+S3 --> C3
##   8.    C3 --> P2+S3
##   9.    C3 --> P2+P3


##           (P,S,C,P,S,C,P,S,C,P)
##            0 1 1 1 2 2 2 3 3 3
S6 <- rbind(c(1,1,0,0,0,0,0,0,0,0),
            c(0,0,1,0,0,0,0,0,0,0),
            c(0,0,1,0,0,0,0,0,0,0),
            c(0,0,0,1,1,0,0,0,0,0),
            c(0,0,0,0,0,1,0,0,0,0),
            c(0,0,0,0,0,1,0,0,0,0),
            c(0,0,0,0,0,0,1,1,0,0),
            c(0,0,0,0,0,0,0,0,1,0),
            c(0,0,0,0,0,0,0,0,1,0))

##           (P,S,C,P,S,C,P,S,C,P)
##            0 1 1 1 2 2 2 3 3 3
P6 <- rbind(c(0,0,1,0,0,0,0,0,0,0),
            c(1,1,0,0,0,0,0,0,0,0),
            c(1,0,0,1,0,0,0,0,0,0),
            c(0,0,0,0,0,1,0,0,0,0),
            c(0,0,0,1,1,0,0,0,0,0),
            c(0,0,0,1,0,0,1,0,0,0),
            c(0,0,0,0,0,0,0,0,1,0),
            c(0,0,0,0,0,0,1,1,0,0),
            c(0,0,0,0,0,0,1,0,0,1))


orderings6 <- findOrderings(source.complexes = S6, product.complexes = P6)



## [7] Example 4.9 (Population dynamics) ###########################################################

## ordered species: 
##   1. A
##   2. B

## ordered reactions:
##   1.   0 --> A
##   2.   A --> 0
##   3.   A --> 2A
##   4. A+B --> B
##   5.   0 --> B
##   6.   B --> 0
##   7.   B --> 2B
##   8. A+B --> A


##           (A,B)
S7 <- rbind(c(0,0),
            c(1,0),
            c(1,0),
            c(1,1),
            c(0,0),
            c(0,1),
            c(0,1),
            c(1,1))

##           (A,B)
P7 <- rbind(c(1,0),
            c(0,0),
            c(2,0),
            c(0,1),
            c(0,1),
            c(0,0),
            c(0,2),
            c(1,0))


orderings7 <- findOrderings(source.complexes = S7, product.complexes = P7)



## [8] Example 4.10 (Histone modification circuit) #################################################

## ordered species: 
##   1. D
##   2. R
##   3. A
##   4. P

## ordered reactions:
##   1.   R --> D
##   2.   D --> R
##   3.   D --> A
##   4.   A --> D
##   5.   A --> A+P
##   6.   P --> 0
##   7. D+P --> A+P
##   8. R+A --> D+R
##   9. D+R --> 2R
##  10. R+A --> D+A
##  11. D+A --> 2A


##           (D,R,A,P)
S8 <- rbind(c(0,1,0,0),
            c(1,0,0,0),
            c(1,0,0,0),
            c(0,0,1,0),
            c(0,0,1,0),
            c(0,0,0,1),
            c(1,0,0,1),
            c(0,1,1,0),
            c(1,1,0,0),
            c(0,1,1,0),
            c(1,0,1,0))

##           (D,R,A,P)
P8 <- rbind(c(1,0,0,0),
            c(0,1,0,0),
            c(0,0,1,0),
            c(1,0,0,0),
            c(0,0,1,1),
            c(0,0,0,0),
            c(0,0,1,1),
            c(1,1,0,0),
            c(0,2,0,0),
            c(1,0,1,0),
            c(0,0,2,0))


orderings8 <- findOrderings(source.complexes = S8, product.complexes = P8)



## [9] Example 4.11 (Ergodicity detection) #########################################################

## ordered species: 
##   1. A
##   2. B

## ordered reactions:
##   1.   0 --> A
##   2.   A --> B
##   3.   B --> 0
##   4. A+B --> A


##           (A,B)
S9 <- rbind(c(0,0),
            c(1,0),
            c(0,1),
            c(1,1))

##           (A,B)
P9 <- rbind(c(1,0),
            c(0,1),
            c(0,0),
            c(1,0))


orderings9 <- findOrderings(source.complexes = S9, product.complexes = P9)


