fn main() {
    println!("Hello, world! Euh, Thomas !");

    // array declaration basics
    const n:usize = 10; //columns
    const m:usize = 5;  // rows

    // mutable float64 array of size (rows = m, columns = n)
    let mut myarray = [[0f64; n]; m];

    for i in 0..m {
        for j in 0..n {
            println!("i: {}, j: {}, array2D[i][j]: {}", i,j, myarray[i][j] )
        };
    };

    for i in 0..m{
        myarray[i][i] = 1f64;
        myarray[i][i+5] = 1f64;
    };
    for i in 0..m {
        for j in 0..n {
            print!(" {}", myarray[i][j] )
        };
        println!(" ");
    };

}
