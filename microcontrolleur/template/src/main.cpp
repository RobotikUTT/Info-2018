//#include "ros.h"
#include <../include/Interface.h>
#include <../include/Verification.h>
#include <../include/Module.h>
#include <../include/Watcher.h>
#include <unistd.h>
#include <time.h>
#include <iostream>

 //ros::NodeHandle nh;
int main(){

  Interface::CreateSingleton();
  Verification::CreateSingleton();
  Module::CreateSingleton();
  Watcher::CreateSingleton();

  while (1) {
    std::cout << "time: " << time(0) << std::endl;
    g_interface_instance_ptr->Update();
    g_verification_instance_ptr->Update();
    g_module_instance_ptr->Update();
    g_watcher_instance_ptr->Update();
    usleep(1000000);
  }
}
