using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;
using steganography.Methods;
using steganography.Functions;
using static System.Net.Mime.MediaTypeNames;
using System.Drawing.Imaging;

namespace steganography
{

    
    static class Program
    {
        /// <summary>
        /// Главная точка входа для приложения.
        /// </summary>
        [STAThread]
        static void Main()
        {
            Bitmap image = new Bitmap(@"C:\Users\Kirill\source\repos\steganography-master\steganography\lena.bmp");
            StegoBitmap stegoImage = new StegoBitmap(image);
            StegoBitmap encrypted = Benham.Hide(stegoImage, "Hello from Kirill Lavrinenko!", Colours.Red);
            encrypted.SaveBitmap(@"C:\Users\Kirill\source\repos\steganography-master\steganography\output.bmp");
            if (Benham.IsHiddenText(encrypted, Colours.Red))
            {
                string result = Benham.GetHiddenText(encrypted, Colours.Red);
                Console.WriteLine("Decoded message:" + result);
            }

            byte[] arr = stegoImage.GetColour(Colours.Red);
            byte[] new_arr = encrypted.GetColour(Colours.Red);

            double PSNR = CommonFunc.calculatePSNR(arr, new_arr,Colours.Red);

            Console.WriteLine("PSNR = " + PSNR.ToString());

            Benham.GetJPEG();
            Bitmap recovered = Benham.ConvertToBitmap(@"C:\Users\Kirill\source\repos\steganography-master\steganography\outputHundred.jpg");
            recovered.Save(@"C:\Users\Kirill\source\repos\steganography-master\steganography\restored.bmp");
            StegoBitmap recovered_JPEG = new StegoBitmap(recovered);

            List<byte> first = Benham.GetBites(encrypted, Colours.Red);
            List<byte> second = Benham.GetBites(recovered_JPEG, Colours.Red);

           int errors = Benham.CalculateErrors(first, second);
            Console.WriteLine("Errors: " + errors);

        }
    }
}
