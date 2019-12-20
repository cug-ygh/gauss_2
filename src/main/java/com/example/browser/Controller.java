package com.example.browser;


import org.springframework.ui.Model;
import org.springframework.web.bind.annotation.RequestMapping;

@org.springframework.stereotype.Controller
@RequestMapping("gauss")
public class Controller {
    private String L;
    private String B;

    @RequestMapping("")
    public String index(){
        return "index";
    }
    @RequestMapping("/go")
    public String go(String s1 ,String s2,String s3,String s4,Model model){
        double res[]=server.GeotoGauss(Double.parseDouble(s1),Double.parseDouble(s2),
                Integer.parseInt(s3),Integer.parseInt(s4));
        model.addAttribute("x",res[0]);
        model.addAttribute("y",res[1]);
        return "go";
    }
    @RequestMapping("/back")
    public String back(String s1 ,String s2,String s3,String s4,String s5,Model model){
        double res[]=server.Gauss2Geo(Double.parseDouble(s1),Double.parseDouble(s2)
                ,Integer.parseInt(s3),Integer.parseInt(s4),Integer.parseInt(s5));
        model.addAttribute("L",res[0]);
        model.addAttribute("B",res[1]);
        return "back";
    }


}
