//#include "ros.h"
#include "../include/Watcher.h"
#include "stdio.h"
#include <iostream>
#include <string>

//include msg

Watcher* watcher_instance_ptr = 0;

Watcher::Watcher(){
}

Watcher::~Watcher(){
  delete watcher_instance_ptr;
}

void Watcher::CreateSingleton(){
  if (watcher_instance_ptr == 0){
    watcher_instance_ptr = new Watcher();
  }
}

void Watcher::Update(){
  std::cout << "Watcher: " << watcher_instance_ptr << "\n" << std::endl;

}
