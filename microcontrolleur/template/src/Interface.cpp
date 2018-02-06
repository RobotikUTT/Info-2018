//#include "ros.h"
#include "../include/Interface.h"
#include "stdio.h"
#include <iostream>
#include <string>

//include msg


Interface* interface_instance_ptr = 0;

Interface::Interface(){
}

Interface::~Interface(){
  delete interface_instance_ptr;
}

void Interface::CreateSingleton(){
  if (interface_instance_ptr == 0){
    interface_instance_ptr = new Interface();
  }
}

void Interface::Update(){
  std::cout << "Interface: " << interface_instance_ptr  << std::endl;
}
