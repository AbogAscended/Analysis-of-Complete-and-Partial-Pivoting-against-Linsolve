# Statistical Analysis of partial and complete pivoting
In this mini-project I wanted to see how different the error between calculating Ax=b was for complete and partial pivoting. I decided to use matlabs linsolve as the actual solution even though its not perfect and calcualte the relative error between
the self implemented methods and linsolve with $\frac{||mine - linsolve||}{||linsolve||}$
