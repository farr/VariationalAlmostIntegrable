open OUnit
open Newton

let test_sqrt () = 
  for i = 1 to 100 do 
    let x2 = Random.float 1.0 in 
    let eps = 100.0*.epsilon_float in 
    let x = 
      newton_solve 10 eps 0.0 eps (fun x -> x*.x -. x2) (fun x -> 2.0*.x) 1.0 in 
    assert_bool "sqrt^2 inaccurate" (abs_float (x*.x -. x2) < 1e4*.epsilon_float)
  done

let tests = "newton.ml tests" >:::
  ["newton's method sqrt" >:: test_sqrt]
