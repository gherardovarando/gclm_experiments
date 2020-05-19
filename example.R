## Example for the paper "Graphical continuous Lyapunov models"
## by Gherardo Varando and Niels Richard Hansen

d <- 5
B <- matrix(
	c( -1,    1,    0,    0,   0,
		 -1,    0,  0.2,    0,   0,
		  0,    0,   -1,    -0.5,   0,
		  0,    0,    0,   -1,   1,
		  0,    0,    1,    0,  -1),
	d, d, byrow = TRUE)
adj <- B != 0

biadj <- diag(d)
#biadj[4, 5] <- biadj[5, 4] <- 1
I <- diag(1, d)
C <- diag(1, d)

round(eigen(B)$values, 2)


BB <- I %x% B + B %x% I
Sigma <- matrix(solve(BB, - as.vector(C)), d, d)
dd <- 1:(d-1)

- zapsmall(B[dd,dd] %*% Sigma[dd, dd] + Sigma[dd, dd] %*% t(B[dd, dd]))

Ctilde <- B[dd, d, drop = FALSE] %*% Sigma[d, dd, drop = FALSE] + 
	Sigma[dd, d, drop = FALSE] %*% t(B[dd, d, drop = FALSE]) + diag(d-1)

zapsmall(Ctilde)

eigen(Ctilde)$values
eigen(B[dd,dd])$values

round(Ctilde, digits = 2)
