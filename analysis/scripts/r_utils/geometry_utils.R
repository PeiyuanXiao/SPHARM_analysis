# geometry_utils.R
# Rotation helper shared by align_svd.R and align_lin2024.R.

# Rotation matrix that maps unit vector `a` onto unit vector `b` (Rodrigues).
get_rot_matrix <- function(a, b) {
  a <- a / sqrt(sum(a^2))
  b <- b / sqrt(sum(b^2))
  cos_theta <- sum(a * b)
  
  # Antiparallel: rotate 180 deg about any axis perpendicular to `a`.
  if (cos_theta < -1 + 1e-10) {
    perp <- if (abs(a[1]) < 0.9) c(1, 0, 0) else c(0, 1, 0)
    v <- perp - sum(perp * a) * a
    v <- v / sqrt(sum(v^2))
    return(2 * outer(v, v) - diag(3))
  }
  
  # Already aligned.
  if (cos_theta > 1 - 1e-10) return(diag(3))
  
  v <- c(
    a[2] * b[3] - a[3] * b[2],
    a[3] * b[1] - a[1] * b[3],
    a[1] * b[2] - a[2] * b[1]
  )
  v_skew <- matrix(
    c( 0,    -v[3],  v[2],
       v[3],  0,    -v[1],
       -v[2],  v[1],  0   ),
    3, 3, byrow = TRUE
  )
  diag(3) + v_skew + v_skew %*% v_skew * ((1 - cos_theta) / sum(v^2))
}