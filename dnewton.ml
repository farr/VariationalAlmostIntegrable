(** Newton solver using DFloat. *)

open Diff.DFloat

let abs_float x = 
  sqrt (x*x)

(** [newton_solve nmax epsabs epsrel epsf f guess] returns the
    argument [x] which makes [f x] zero.  The function uses Newton
    iteration.  The iteration stops when either

    - [f x < epsf]

    - [abs_float (x_new -. x_old) < epsabs +. epsrel*.(abs_float x_new)] 

    - The number of iterations is greater than [nmax], in which case
      [Failure] is raised. *)
let newton_solve nmax epsabs epsrel epsf f x0 = 
  let f' = d f in 
  let rec loop i x = 
    if Pervasives.(>=) i nmax then 
      raise (Failure "newton_solve: too many iterations")
    else
      let fx = f x and 
          f'x = f' x in 
      if f'x = (const 0.0) then 
        raise (Failure "newton_solve: zero slope encountered")
      else if abs_float fx < (const epsf) then 
        x
      else
        let x_new = x - fx/f'x in 
        if abs_float (x - x_new) < (const epsabs) + (const epsrel)*(abs_float x_new) then 
          x_new
        else
          loop (Pervasives.(+) i 1) x_new in 
  loop 0 x0

let same_sign a b = 
  ((a < zero) && (b < zero)) ||
  ((a > zero) && (b > zero))

let bracket_x xmax xmin fxmax fxmin = 
  (fxmax*xmin - fxmin*xmax)/(fxmax - fxmin)

let between xmin x xmax = 
  (xmin <= x) && (xmax >= x)

let bracket_solve nmax epsabs epsrel epsf f xmin xmax = 
  let df = d f in 
  let rec loop i xmin xmax fxmin fxmax x0 = 
    if Pervasives.(>) i nmax then 
      raise (Failure "bracket_solve: too many iterations")
    else
      let fx0 = f x0 and 
          dfx0 = df x0 in 
      let xn = x0 - fx0/dfx0 in 
      let x1 = if between xmin xn xmax then xn else (const 0.5)*(xmin+xmax) in 
      let fx1 = f x1 in 
      let dx = abs_float (x0 - x1) and 
          xmag = abs_float x1 in 
      if (abs_float fx1 < (const epsf)) || (dx < (const epsabs) + (const epsrel)*xmag) then 
        x1
      else
        let (xmin,fxmin,xmax,fxmax) = 
          if same_sign fx1 fxmin then 
            (x1,fx1,xmax,fxmax)
          else
            (xmin,fxmin,x1,fx1) in 
        loop (Pervasives.(+) i 1) xmin xmax fxmin fxmax x1 in 
  loop 0 xmin xmax (f xmin) (f xmax) ((const 0.5)*(xmin+xmax))
