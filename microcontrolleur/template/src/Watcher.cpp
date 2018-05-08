//#include "ros.h"
#include "../include/Watcher.h"
#include "stdio.h"
#include <iostream>
#include <string>

//include msg

Watcher* g_watcher_instance_ptr = 0;

Watcher::Watcher(){
}

Watcher::~Watcher(){
  delete g_watcher_instance_ptr;
}

void Watcher::CreateSingleton(){
  if (g_watcher_instance_ptr == 0){
    g_watcher_instance_ptr = new Watcher();
  }
}

void Watcher::Update(){
  std::cout << "Watcher: " << g_watcher_instance_ptr << "\n" << std::endl;

}
