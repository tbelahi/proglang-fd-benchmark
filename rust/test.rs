fn main(){

    const N:usize = 1000;
    const M:usize = 1000;

    //let mut array = [[5f64; N]; M];

    //println!("{}", array[0][0]);

    let mut _2d = Vec::new();
    _2d.push(Vec::new());
    match _2d.last_mut() {
        Some(v) => v.push(1f64),
        None => ()
    }
    println!("{}", _2d[0][0]);

    let mut row = std::iter::repeat(0f64).take(N).collect();
    let mut pseudo_array: Vec<Vec<f64>> = std::iter::repeat(row).take(M).collect();

    println!("constructed");
    println!("pseudo array, en bas à droite {:?}", pseudo_array[M-1][N-1]);

}
