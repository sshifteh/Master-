_new[i][j][k] = u[i][j][ĸ]
+ dt*(nu*u[i-1][j][k]+ nu*u[i][j][k-1] + nu*u[i][j-1][k] - 6*nu*u[i][j][k]
+ nu*u[i][j+1][k] + nu*u[i][j][k+1]+ nu*u[i+1][j][k] + f[i][j][k])
