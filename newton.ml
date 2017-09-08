(** Newton solver. *)

(** [newton_solve nmax epsabs epsrel epsf f f' guess] returns the
    argument [x] which makes [f x] zero.  The function uses Newton
    iteration, and therefore requires knowledge of [f'], the
    derivative of [f].  The iteration stops when either

    - [f x < epsf]

    - [abs_float (x_new -. x_old) < epsabs +. epsrel*.(abs_float x_new)] 

    - The number of iterations is greater than [nmax], in which case
      [Failure] is raised. *)
let newton_solve nmax epsabs epsrel epsf f f' x0 = 
  let rec loop i x = 
    if i >= nmax then 
      raise (Failure "newton_solve: too many iterations")
    else
      let fx = f x and 
          f'x = f' x in 
      if f'x = 0.0 then 
        raise (Failure "newton_solve: zero slope encountered")
      else if abs_float fx < epsf then 
        x
      else
        let x_new = x -. fx/.f'x in 
        if abs_float (x -. x_new) < epsabs +. epsrel*.(abs_float x_new) then 
          x_new
        else
          loop (i+1) x_new in 
  loop 0 x0
