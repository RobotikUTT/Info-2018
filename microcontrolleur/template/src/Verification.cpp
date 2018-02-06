//#include "ros.h"
#include "../include/Verification.h"
#include "stdio.h"
#include <iostream>
#include <string>

//include msg

Verification* verification_instance_ptr = 0;


Verification::Verification(){
}

Verification::~Verification(){
  delete verification_instance_ptr;
}

void Verification::CreateSingleton(){
  if (verification_instance_ptr == 0){
    verification_instance_ptr = new Verification();
  }
}

void Verification::Update(){
  std::cout << "Verification: " << verification_instance_ptr << std::endl;
}
