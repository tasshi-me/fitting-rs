fn main() {
    println!("Hello, world!");
}

#[test]
fn it_works() {}

#[test]
#[should_panic]
fn it_does_not_works() {
    assert!(false);
}
