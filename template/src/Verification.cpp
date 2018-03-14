//#include "ros.h"
#include "../include/Verification.h"
#include "stdio.h"
#include <iostream>
#include <string>

//include msg

Verification* g_verification_instance_ptr = 0;


Verification::Verification(){
}

Verification::~Verification(){
  delete g_verification_instance_ptr;
}

void Verification::CreateSingleton(){
  if (g_verification_instance_ptr == 0){
    g_verification_instance_ptr = new Verification();
  }
}

void Verification::Update(){
  std::cout << "Verification: " << g_verification_instance_ptr << std::endl;
}
