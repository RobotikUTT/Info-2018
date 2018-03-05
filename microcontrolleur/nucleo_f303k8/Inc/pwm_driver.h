#ifndef __PWM_DRIVER_H__
#define __PWM_DRIVER_H__

#include "stm32f3xx_hal.h"
//#include <stdio.h>
//#include "stm32f3xx_hal_dma.h"
//#include "stm32f3xx_hal_tim.h"

void set_timer_freq(TIM_HandleTypeDef* timer,uint32_t freq);
void set_channel_duty_cycle(TIM_HandleTypeDef* timer, 
                uint32_t Channel, uint8_t duty_cycle);


#endif