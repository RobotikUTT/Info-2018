//#include "ros.h"
#include "../include/Module.h"
#include "stdio.h"
#include <iostream>
#include <string>

//include msg

Module* module_instance_ptr = 0;

Module::Module(){
}

Module::~Module(){
  delete module_instance_ptr;
}

void Module::CreateSingleton(){
  if (module_instance_ptr == 0){
    module_instance_ptr = new Module();
  }
}

void Module::Update(){
  std::cout << "Module: " << module_instance_ptr << std::endl;
  verification_instance_ptr->Update();

}
