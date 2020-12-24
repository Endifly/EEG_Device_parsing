*** ไฟล์ _gains ***
ใช้สำหรับเก็บข้อมูล gain ของแต่ละ device โดยแต่ละไฟล์จะมี format

[device_name] {
	"ch" : {
		"AF7" : [([freq],[gain], ...] ,
		"AF8" : [([freq],[gain], ...] ,
	} ,
	"params" : {
      "amplitude" : {
        "value" : [profile.func_gen_amp_l],
        "unit" : "V"
      },
      "offset" : {
        "value" :  [profile.func_gen_offs_l],
        "unit" : "V"
      },
      "duty" : {
        "value" : [profile.func_gen_duty_l],
        "unit" : None
      },
      "phase" : {
        "value" : [profile.func_gen_phase_l],
        "unit" : "degree"
      }, 
      "desc" : "Gain of amplitude "
    }

*** Device Data ***	
ไฟล์ข้อมูล Raw Data ของแต่ละ Device ให้ตั้งชื่อไฟล์เริ่มต้นด้วยตัวพิมพ์ใหญ่ และภายในเก็บ Raw Data เอาไว้
- หากไม่มี Folder ย่อย profile_gain_compute จะให้ชื่อ Channel ว่า AF7
- หากมี Folder ย่อย profile_gain_compute จะใช้ชื่อ Folder ที่เก็บไฟล์ Raw Data เป็นชื่อ Chaneel


***	profile_loader ***	
ใช้สำหรับอ่านข้อมูล Device Data และ แสดง FFT_PW, FFT_AMP ในแต่ละ Freq

***	profile_gain_compute ***	
ใช้สำหรับอ่าน Device Data และคำนวน gain ในแต่ละ Freq และเก็บข้อมูลลง _gains

***	profile_device_parser ***	
ใช้สำหรับอ่าน Device Data และ _gains และแปลงข้อมูลลักษณะ Amplitude จาก master_device ไปหา device อื่น
