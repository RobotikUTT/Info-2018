//#include "ros.h"
#include "../include/Interface.h"
#include "stdio.h"
#include <iostream>
#include <string>

//include msg


Interface* g_interface_instance_ptr = 0;

Interface::Interface(){
}

Interface::~Interface(){
  delete g_interface_instance_ptr;
}

void Interface::CreateSingleton(){
  if (g_interface_instance_ptr == 0){
    g_interface_instance_ptr = new Interface();
  }
}

void Interface::Update(){
  std::cout << "Interface: " << g_interface_instance_ptr  << std::endl;
}
