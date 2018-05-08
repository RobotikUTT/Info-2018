//#include "ros.h"
#include "../include/Module.h"
#include "stdio.h"
#include <iostream>
#include <string>

//include msg

Module* g_module_instance_ptr = 0;

Module::Module(){
}

Module::~Module(){
  delete g_module_instance_ptr;
}

void Module::CreateSingleton(){
  if (g_module_instance_ptr == 0){
    g_module_instance_ptr = new Module();
  }
}

void Module::Update(){
  std::cout << "Module: " << g_module_instance_ptr << std::endl;
  g_verification_instance_ptr->Update();

}
