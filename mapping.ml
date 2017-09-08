(** Procedures to implement mapping-type integrators. *)

(** [fold_righti f arr init] returns [f 0 arr.(0) (f 1 arr.(1) (... (f
    (n-1) arr.(n-1) init)))]*)
let fold_righti f arr init = 
  let state = ref init in 
  for i = Array.length arr - 1 downto 0 do 
    state := f i arr.(i) !state
  done;
  !state

(** [mapping_advance cs a b dt state] returns the composition of [a
    (cs.(0)*.dt)], [b (cs.(1)*.dt)], ..., acting on [state], where [a]
    and [b] are evolution operators for the [state] *)
let mapping_advance coeffs eA eB dt state = 
  fold_righti 
    (fun i coeff state -> 
      if i mod 2 = 0 then 
        eA (dt *. coeff) state
      else
        eB (dt *. coeff) state)
    coeffs
    state
