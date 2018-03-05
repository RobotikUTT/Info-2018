#include "pwm_driver.h"

void set_timer_freq(TIM_HandleTypeDef* timer,uint32_t freq)
{
  uint32_t sys_freq = HAL_RCC_GetPCLK1Freq();

  timer->Init.Prescaler = 0;
  timer->Init.CounterMode = TIM_COUNTERMODE_UP;
  timer->Init.Period = (uint16_t)(2*sys_freq/(freq));
  timer->Init.ClockDivision = TIM_CLOCKDIVISION_DIV1;
  timer->Init.AutoReloadPreload = TIM_AUTORELOAD_PRELOAD_DISABLE;
  if (HAL_TIM_PWM_Init(timer) != HAL_OK)
  {
    _Error_Handler(__FILE__, __LINE__);
  }
}

void set_channel_duty_cycle(TIM_HandleTypeDef* timer, 
                uint32_t Channel, uint8_t duty_cycle)
{
	TIM_OC_InitTypeDef sConfigOC;
	sConfigOC.OCMode = TIM_OCMODE_PWM1;
	sConfigOC.Pulse = timer->Init.Period*duty_cycle/0xFF;
	sConfigOC.OCPolarity = TIM_OCPOLARITY_HIGH;
	sConfigOC.OCFastMode = TIM_OCFAST_DISABLE;
	if (HAL_TIM_PWM_ConfigChannel(timer, &sConfigOC, Channel) != HAL_OK)
	{
	  _Error_Handler(__FILE__, __LINE__);
	}
}