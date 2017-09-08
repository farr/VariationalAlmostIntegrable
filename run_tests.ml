open OUnit

let _ = Randomize.init ()

let tests = "all tests" >:::
  ["sho.ml tests" >: Sho_test.tests;
   "newton.ml tests" >: Newton_test.tests;
   "dthree_body.ml tests" >: Dthree_body_test.tests]

let _ = 
  let results = run_test_tt_main tests in 
  let nfail = 
    List.fold_left
      (fun n result -> 
        match result with 
        | RSuccess(_) -> n
        | _ -> n+1)
      0
      results in 
  exit nfail
