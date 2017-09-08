open OUnit
open Diff.DFloat
open Dthree_body

let random_vector () = 
  Array.init 3 (fun _ -> C (Random.float 1.0 -. 0.5))

let test_jacobi_lagrangian () = 
  for i = 1 to 100 do 
    let sys = 
      ({m = C (Random.float 1.0);
        q = random_vector ();
        v = random_vector ()}, 
       {m = C (Random.float 1.0);
        q = random_vector ();
        v = random_vector ()},
       {m = C (Random.float 1.0);
        q = random_vector ();
        v = random_vector ()}) in 
    let l = 
      match lagrangian sys with 
      | C x -> x
      | _ -> raise (Failure "lagrangian didn't produce a number") and 
        lj = 
      match jacobi_lagrangian (jacobi_transform sys) with 
      | C x -> x
      | _ -> raise (Failure "jacobi_lagrangian didn't produce a numebr") in 
    assert_bool "l and lj not equal" (Pervasives.(<) (abs_float (l -. lj)) 1e-8)
  done

let tests = "dthree_body.ml tests" >:::
  ["jacobi_lagrangian tests" >:: test_jacobi_lagrangian]
