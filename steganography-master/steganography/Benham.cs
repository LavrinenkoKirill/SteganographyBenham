using System;
using System.Collections.Generic;
using System.Diagnostics.Eventing.Reader;
using System.Drawing;
using System.Drawing.Imaging;
using System.IO;
using System.Text;
using steganography.Functions;

namespace steganography.Methods
{
    class Benham
    {
        private const int SizeSegment = 8; // размер сегмента
        private const int CoefDif = 50; // разность коэффициентов
        public const double MinHighFreq = 40; // порог для высокочастотных коэффициентов
        public const double MaxLowFreq = 2600;  // порог для низкочастотных коэффициентов

        static Point p1 = new Point(6, 3); // точки для работы с сегментами коэффициентов ДКП 
        static Point p2 = new Point(4, 4);
        static Point p3 = new Point(3, 6); 

        
        public static StegoBitmap Hide(StegoBitmap stgBitmap, string txt, Colours colour)
        {
            Bitmap bmap = stgBitmap.GetImage();
            var arrForHiding = new byte[bmap.Width, bmap.Height];
            if ((bmap.Width % SizeSegment) != 0 || (bmap.Height % SizeSegment) != 0) //если изображение не делится на сегменты ровно 
                CommonFunc.Cut(ref bmap, SizeSegment);
            for (int i = 0; i < bmap.Width; i++)
            {
                for (int j = 0; j < bmap.Height; j++)
                {
                    if (colour == Colours.Red)
                        arrForHiding[i, j] = bmap.GetPixel(i, j).R;
                    else if (colour == Colours.Green)
                        arrForHiding[i, j] = bmap.GetPixel(i, j).G;
                    else if (colour == Colours.Blue)
                        arrForHiding[i, j] = bmap.GetPixel(i, j).B;
                    else
                        throw new NullReferenceException();
                }
            }
            var stgByte = Encoding.GetEncoding(1251).GetBytes(txt); 
            byte[] len = CommonFunc.LenInBytes(txt.Length, CommonFunc.Size(txt.Length));       
            byte[] txtByte = new byte[2 + len.Length + stgByte.Length];
            txtByte[0] = Convert.ToByte('B'); 
            txtByte[1] = CommonFunc.Size(txt.Length); 
            int ind = 0;
            for (int i = 0; i < len.Length; i++)
                txtByte[i + 2] = len[ind++];  
            ind = 0;
            for (int i = 0; i < stgByte.Length; i++)
                txtByte[i + 2 + len.Length] = stgByte[ind++];
            var segm = new List<byte[,]>();
            Separate(arrForHiding, segm, bmap.Width, bmap.Height, SizeSegment);// разбиваем массив на сегменты
            // дискретное косинусное преобразование
            var dctList = new List<double[,]>();
            foreach (var b in segm)
                dctList.Add(DCT(b)); // список из коэффициентов ДКП
            
            SetText(txtByte, ref dctList, CoefDif); // внедрение текста
            // обратное дискретное косинусное преобразование
            var idctList = new List<double[,]>();
            foreach (var d in dctList)
                idctList.Add(IDCT(d));
            var newArr = new double[bmap.Width, bmap.Height]; 
            Join(ref newArr, idctList, bmap.Width, bmap.Height, SizeSegment); //соединяем сегменты
            var normArr = Normalize(newArr); // модификация коэффициентов ДКП иногда приводит к выходу значений интенсивностей пикселей изображения за пределы допустимого диапазона [0,255], проводим нормирование указанных значений
            return new StegoBitmap(bmap, normArr, colour);
        }

        // разбиение массива на сегменты
        private static void Separate(byte[,] arr, List<byte[,]> segmList, int width, int height, int sizeSegm)
        {
            int numSW = width / sizeSegm; 
            int numSH = height / sizeSegm;
            for (int i = 0; i < numSW; i++)
            {
                int firstWPoint = i * sizeSegm; 
                int lastWPoint = firstWPoint + sizeSegm - 1; 
                for (int j = 0; j < numSH; j++)
                {
                    int firstHPoint = j * sizeSegm; 
                    int lastHPoint = firstHPoint + sizeSegm - 1; 
                    segmList.Add(SegmBytes(arr, firstWPoint, lastWPoint, firstHPoint, lastHPoint));
                }
            }
        }

        
        private static byte[,] SegmBytes(byte[,] arr, int a, int b, int c, int d)
        {
            var sg = new byte[b - a + 1, d - c + 1];
            for (int i = a, x = 0; i <= b; i++, x++)
                for (int j = c, y = 0; j <= d; j++, y++)
                    sg[x, y] = arr[i, j];
            return sg;
        }

        // определение значений коэффициентов для текущего значения аргументов
        private static double GetCoefficient(int arg)
        {
            if (arg == 0)
                return 1.0 / Math.Sqrt(2);
            return 1;
        }

        private static double[,] DCT(byte[,] b)
        {
            int len = b.GetLength(0); // получаем размер сегмента
            double[,] arrDCT = new double[len, len]; // новый массив после дискретного косинусного преобразования
            double temp = 0;
            for (int i = 0; i < len; i++)
            {
                for (int j = 0; j < len; j++)
                {
                    temp = 0;
                    for (int x = 0; x < len; x++)
                    {
                        for (int y = 0; y < len; y++)
                            temp += b[x, y] * Math.Cos(Math.PI * i * (2 * x + 1) / (2 * len)) * Math.Cos(Math.PI * j * (2 * y + 1) / (2 * len)); // вычисление спектральных коэффициентов ДКП для каждого сегмента (без умножения на значения коэффициентов для текущего значения аргумента) 
                    }
                    arrDCT[i, j] = GetCoefficient(j) * GetCoefficient(i) * temp / Math.Sqrt(2 * len); // вычисление спектральных коэффициентов ДКП для каждого сегмента (домножаем на значения коэффициентов для текущего значения аргумента)
                }
            }
            return arrDCT;
        }
        private static double[,] DCT(double [,] b)
        {
            int len = b.GetLength(0); // получаем размер сегмента
            double[,] arrDCT = new double[len, len]; // новый массив после дискретного косинусного преобразования
            double temp = 0;
            for (int i = 0; i < len; i++)
            {
                for (int j = 0; j < len; j++)
                {
                    temp = 0;
                    for (int x = 0; x < len; x++)
                    {
                        for (int y = 0; y < len; y++)
                            temp += b[x, y] * Math.Cos(Math.PI * i * (2 * x + 1) / (2 * len)) * Math.Cos(Math.PI * j * (2 * y + 1) / (2 * len)); // вычисление спектральных коэффициентов ДКП для каждого сегмента (без умножения на значения коэффициентов для текущего значения аргумента) 
                    }
                    arrDCT[i, j] = GetCoefficient(j) * GetCoefficient(i) * temp / Math.Sqrt(2 * len); // вычисление спектральных коэффициентов ДКП для каждого сегмента (домножаем на значения коэффициентов для текущего значения аргумента)
                }
            }
            return arrDCT;
        }

        //Блоки, не отвечающие требованию 1, характеризуются наличием слишком больших значений низкочастотных коэффициентов ДКП, сопоставимых по своей величине с DC-коэффнциентом.
        // проверка резкости переходов яркости
        public static bool ValidSharpness(double[,] dl)
        {
            double lowFreq = 0;
            int segmEnd = SizeSegment - 1;
            for (int i = 0; i < SizeSegment - 1; i++)
            {
                for (int j = 0; j < segmEnd; j++)
                    lowFreq += Math.Abs(dl[i, j]); 
                segmEnd--;
            }
            lowFreq -= Math.Abs(dl[0, 0]); 
            if (MaxLowFreq - lowFreq >= 0.001)
                return true;
            else
                return false;
        }

        //Для блоков, не удовлетворяющих требованию 2, характерно равенство нулю большинства высокочастотных коэффициентов. 
        // проверка монотонности
        public static bool ValidMoonotony(double [,] dl)
        {
            double highFreq = 0;
            int segmEnd = SizeSegment;
            for (int i = 2; i < SizeSegment; i++)
            {
                for (int j = segmEnd - 1; j < SizeSegment; j++)
                    highFreq += Math.Abs(dl[i, j]);
                segmEnd--;
            }
            if (highFreq - MinHighFreq >= 0.001)
                return true;
            else
                return false;
        }

        //проверка, что после обратного ДКП и повторного ДКП блок также остается пригодным для внедрения
        public static bool TestValidInver(double[,] dl)
        {
            var idct_dl = IDCT(dl);            
            var byte_dl  = Normalize(idct_dl);
            var dct_dl = DCT(byte_dl);
            if (ValidMoonotony(dct_dl) && ValidSharpness(dct_dl))
                return true;
            else
                return false;
        }

        
        // внедрение текста
        private static void SetText(byte[] txt, ref List<double[,]> blocks, int coefDif)
        {
            List<int> freePos = new List<int>(); // свободные позиции в соответствии с размером списка коэффициентов ДКП
            for (int i = 0; i < blocks.Count; i++)
                freePos.Add(i);
            for (int i = 0; i < txt.Length; i++)
            {
                bool[] bitsSymb = CommonFunc.ByteBoolArr(txt[i]); 
                for (int j = 0; j < 8; j++)
                {
                    bool currentBit = bitsSymb[j];
                    int pos = freePos[0]; 
                    freePos.RemoveAt(0);
                    
                    // берем значения коэффициентов ДКП в точках
                    double P1 = blocks[pos][p1.X, p1.Y];
                    double P2 = blocks[pos][p2.X, p2.Y];
                    double P3 = blocks[pos][p3.X, p3.Y];
                    bool HideBit = false;
                    if (ValidMoonotony(blocks[pos]) && ValidSharpness(blocks[pos]))
                    {
                        if (currentBit) //для передачи бита "1" стремятся, чтобы третий коэффициент стал большим по сравнению с первым и вторым коэффициентами
                        {
                            
                            if ((Math.Max(P1, P2) - P3) > 0.001 || Math.Abs(P3 - Math.Max(P1, P2)) < 0.001)
                            {
                                P3 = Math.Max(P1, P2) + coefDif / 2;
                                if (P1 > P2)
                                    P1 -= coefDif / 2;
                                else                               
                                    P2 -= coefDif / 2;
                            }
                        }
                        else //для передачи бита "0" стремятся, чтобы третий коэффициент стал меньше любого из первых двух
                        {
                            
                            if ((P3 - Math.Min(P1, P2)) > 0.001 || Math.Abs(P3 - Math.Min(P1, P2)) < 0.001)
                            {
                                P3 = Math.Min(P1, P2) - coefDif / 2;
                                if (P1 < P2)
                                    P1 += coefDif / 2;
                                else
                                    P2 += coefDif / 2;
                                
                            }
                        }
                        HideBit = true;
                    }
                    //присваиваем коэффициентам ДКП новые значения
                    blocks[pos][p1.X, p1.Y] = P1;
                    blocks[pos][p2.X, p2.Y] = P2;
                    blocks[pos][p3.X, p3.Y] = P3;
                    if (HideBit && TestValidInver(blocks[pos]))
                        j++;
                    
                    j--;
                }
            }
        }

        // обратное дискретное косинусное преобразование
        private static double[,] IDCT(double[,] dct)
        {
            int len = dct.GetLength(0); // получаем размер сегмента ДКП
            double[,] result = new double[len, len]; // новый массив после обратного дискретного косинусного преобразования
            double temp = 0;
            for (int i = 0; i < len; i++)
            {
                for (int j = 0; j < len; j++)
                {
                    temp = 0;
                    for (int x = 0; x < len; x++)
                    {
                        for (int y = 0; y < len; y++)
                            temp += GetCoefficient(x) * GetCoefficient(y) * dct[x, y] * Math.Cos(Math.PI * x * (2 * i + 1) / (2 * len)) * Math.Cos(Math.PI * y * (2 * j + 1) / (2 * len));
                    }
                    result[i, j] = temp / (Math.Sqrt(2 * len));
                }
            }
            return result;
        }

        // соединяем сегменты
        private static void Join(ref double[,] arr, List<double[,]> Idct, int width, int height, int sizeSegm)
        {
            var temp = Idct.ToArray();
            int numSW = width / sizeSegm; 
            int numSH = height / sizeSegm; 
            int k = 0;
            for (int i = 0; i < numSW; i++)
            {
                int firstWPoint = i * sizeSegm; 
                int lastWPoint = firstWPoint + sizeSegm - 1; 
                for (int j = 0; j < numSH; j++)
                {
                    int firstHPoint = j * sizeSegm;
                    int lastHPoint = firstHPoint + sizeSegm - 1;
                    Insert(ref arr, temp[k], firstWPoint, lastWPoint, firstHPoint, lastHPoint);
                    k++;
                }
            }
        }

        // вставка сегмент в массив
        private static void Insert(ref double[,] arr, double[,] temp, int firstWPoint, int lastWPoint, int firstHPoint, int lastHPoint)
        {
            for (int i = firstWPoint, u = 0; i < lastWPoint + 1; i++, u++)
            {
                for (int j = firstHPoint, v = 0; j < lastHPoint + 1; j++, v++)
                    arr[i, j] = temp[u, v];
            }
        }

        // нормировка
        private static byte[,] Normalize(double[,] idct_dl)
        {
           
            var byte_dl = new byte[idct_dl.GetLength(0), idct_dl.GetLength(1)];
            for (int i = 0; i < idct_dl.GetLength(0); i++)
            {
                for (int j = 0; j < idct_dl.GetLength(1); j++)
                {
                    if (idct_dl[i, j] > 255)
                        byte_dl[i, j] = 255;
                    else if (idct_dl[i, j] < 0)
                        byte_dl[i, j] = 0;
                    else
                        byte_dl[i, j] = Convert.ToByte(idct_dl[i, j]);
                }
            }
            return byte_dl;
        }

        // проверка на налчие скрытого текста
        public static bool IsHiddenText(StegoBitmap stgbmap, Colours c)
        {
            var bmap = stgbmap.GetImage();
            int width = bmap.Width;
            int height = bmap.Height;
            var arrWhereHide = new byte[bmap.Width, bmap.Height]; 
            for (int i = 0; i < bmap.Width; i++)
            {
                for (int j = 0; j < bmap.Height; j++)
                {
                    if (c == Colours.Red)
                        arrWhereHide[i, j] = bmap.GetPixel(i, j).R;
                    else if (c == Colours.Green)
                        arrWhereHide[i, j] = bmap.GetPixel(i, j).G;
                    else if (c == Colours.Blue)
                        arrWhereHide[i, j] = bmap.GetPixel(i, j).B;
                    else
                        throw new NullReferenceException();
                }
            }
            int numSegm = bmap.Width * bmap.Height / (SizeSegment * SizeSegment); //общее число сегментов
            var segm = new List<byte[,]>();
            Separate(arrWhereHide, segm, bmap.Width, bmap.Height, SizeSegment);// разбиваем массив на сегменты
            // дискретное косинусное преобразование
            var dctList = new List<double[,]>();
            foreach (var b in segm)
                dctList.Add(DCT(b)); // список из коэффициентов ДКП
            var txtByte = new List<byte>();
            List<int> possibPos = new List<int>(); // возможные позиции в соответствии с размером списка коэффициентов ДКП
            for (int i = 0; i < dctList.Count; i++)
                possibPos.Add(i);
            var bits = new bool[8]; 
            for (int j = 0; j < 8; j++)
            {
                int pos = possibPos[0]; 
                possibPos.RemoveAt(0);
                if (ValidMoonotony(dctList[pos]) && ValidSharpness(dctList[pos]))
                {
                    double P1 = dctList[pos][p1.X, p1.Y];
                    double P2 = dctList[pos][p2.X, p2.Y];
                    double P3 = dctList[pos][p3.X, p3.Y];
                    
                    if (Math.Min(P1, P2) - P3 >= 0.001)
                        bits[j] = false;
                    else if (P3 - Math.Max(P1, P2) >= 0.001)
                        bits[j] = true;

                    j++;
                }
                j--;
            }
            txtByte.Add(CommonFunc.BoolArrByte(bits));
            if (txtByte.ToArray()[0] != Convert.ToByte('B')) 
                return false;
            return true;
        }

        public static string GetHiddenText(StegoBitmap stgBitmap, Colours c)
        {
            var bmap = stgBitmap.GetImage();
            int width = bmap.Width;
            int height = bmap.Height;
            var arrWhereHide = new byte[bmap.Width, bmap.Height]; 
            for (int i = 0; i < bmap.Width; i++)
            {
                for (int j = 0; j < bmap.Height; j++)
                {
                    if (c == Colours.Red)
                        arrWhereHide[i, j] = bmap.GetPixel(i, j).R;
                    else if (c == Colours.Green)
                        arrWhereHide[i, j] = bmap.GetPixel(i, j).G;
                    else if (c == Colours.Blue)
                        arrWhereHide[i, j] = bmap.GetPixel(i, j).B;
                    else
                        throw new NullReferenceException();
                }
            }
            var segm = new List<byte[,]>();
            Separate(arrWhereHide, segm, bmap.Width, bmap.Height, SizeSegment);// разбиваем массив на сегменты
            // дискретное косинусное преобразование
            var dctList = new List<double[,]>();
            foreach (var b in segm)
                dctList.Add(DCT(b)); // список из коэффициентов ДКП
            var txtByte = new List<byte>();
            
            List<int> possibPos = new List<int>(); // возможные позиции в соответствии с размером списка коэффициентов ДКП
            for (int i = 0; i < dctList.Count; i++)
                possibPos.Add(i);
            int end = 2; //конец прохода
            bool LenDone = true;
                for (int i = 0; i < end; i++)
                {
                    var bits = new bool[8]; // булевый массив битов символа
                    for (int j = 0; j < 8; j++)
                    {
                        int pos = possibPos[0]; // позиция 
                        possibPos.RemoveAt(0);
                        if (ValidMoonotony(dctList[pos]) && ValidSharpness(dctList[pos]))
                        {
                            
                            double P1 = dctList[pos][p1.X, p1.Y];
                            double P2 = dctList[pos][p2.X, p2.Y];
                            double P3 = dctList[pos][p3.X, p3.Y];
                            
                            if (Math.Min(P1, P2) - P3 >= 0.001)
                                bits[j] = false;
                            else if (P3 - Math.Max(P1, P2) >= 0.001)
                                bits[j] = true;

                            j++;
                        }
                        j--;
                    }
                    txtByte.Add(CommonFunc.BoolArrByte(bits));
                    if (i == 0 && txtByte.ToArray()[0] != Convert.ToByte('B')) // в случае отсутвия метки возвращаем пустую строку
                        return "";
                    else if (i == 1) 
                    {
                        end = txtByte.ToArray()[1] + 2;
                        txtByte.Clear();
                    }
                    else if (LenDone && i + 1 == end) 
                    {
                        end += CommonFunc.IntBytes(txtByte);
                        txtByte.Clear();
                        LenDone = false;
                    }
                }
                return Encoding.GetEncoding(1251).GetString(txtByte.ToArray());  
        }


        public static ImageCodecInfo GetEncoder(ImageFormat format)
        {
            ImageCodecInfo[] codecs = ImageCodecInfo.GetImageEncoders();
            foreach (ImageCodecInfo codec in codecs)
            {
                if (codec.FormatID == format.Guid)
                {
                    return codec;
                }
            }
            return null;
        }


        public static void GetJPEG()
        { 
            using (Bitmap bmp1 = new Bitmap(@"C:\Users\Kirill\source\repos\steganography-master\steganography\output.bmp"))
            {
                ImageCodecInfo jpgEncoder = GetEncoder(ImageFormat.Jpeg);
                System.Drawing.Imaging.Encoder myEncoder =
                    System.Drawing.Imaging.Encoder.Quality;
 
                EncoderParameters myEncoderParameters = new EncoderParameters(1);

                EncoderParameter myEncoderParameter = new EncoderParameter(myEncoder, 100L);
          

                myEncoderParameters.Param[0] = myEncoderParameter;
                bmp1.Save(@"C:\Users\Kirill\source\repos\steganography-master\steganography\outputHundred.jpg", jpgEncoder, myEncoderParameters);
            }
        }


        public static Bitmap ConvertToBitmap(string fileName)
        {
            Bitmap bitmap;
            using (Stream bmpStream = System.IO.File.Open(fileName, System.IO.FileMode.Open))
            {
                Image image = Image.FromStream(bmpStream);

                bitmap = new Bitmap(image);

            }
            return bitmap;
        }

        public static List<byte> GetBites(StegoBitmap stgBitmap, Colours c)
        {
            var bmap = stgBitmap.GetImage();
            int width = bmap.Width;
            int height = bmap.Height;
            var arrWhereHide = new byte[bmap.Width, bmap.Height]; 
            for (int i = 0; i < bmap.Width; i++)
            {
                for (int j = 0; j < bmap.Height; j++)
                {
                    if (c == Colours.Red)
                        arrWhereHide[i, j] = bmap.GetPixel(i, j).R;
                    else if (c == Colours.Green)
                        arrWhereHide[i, j] = bmap.GetPixel(i, j).G;
                    else if (c == Colours.Blue)
                        arrWhereHide[i, j] = bmap.GetPixel(i, j).B;
                    else
                        throw new NullReferenceException();
                }
            }
            var segm = new List<byte[,]>();
            Separate(arrWhereHide, segm, bmap.Width, bmap.Height, SizeSegment);
            var dctList = new List<double[,]>();
            foreach (var b in segm)
                dctList.Add(DCT(b)); 
            var txtByte = new List<byte>();

            List<int> possibPos = new List<int>(); 
            for (int i = 0; i < dctList.Count; i++)
                possibPos.Add(i);
            int end = 2; 
            bool LenDone = true;
            for (int i = 0; i < end; i++)
            {
                var bits = new bool[8]; 
                for (int j = 0; j < 8; j++)
                {
                    int pos = possibPos[0];  
                    possibPos.RemoveAt(0);
                    if (ValidMoonotony(dctList[pos]) && ValidSharpness(dctList[pos]))
                    {

                        double P1 = dctList[pos][p1.X, p1.Y];
                        double P2 = dctList[pos][p2.X, p2.Y];
                        double P3 = dctList[pos][p3.X, p3.Y];
                  
                        if (Math.Min(P1, P2) - P3 >= 0.001)
                            bits[j] = false;
                        else if (P3 - Math.Max(P1, P2) >= 0.001)
                            bits[j] = true;

                        j++;
                    }
                    j--;
                }
                txtByte.Add(CommonFunc.BoolArrByte(bits));
                if (i == 1) 
                {
                    end = txtByte.ToArray()[1] + 2;
                    txtByte.Clear();
                }
                else if (LenDone && i + 1 == end) 
                {
                    end += CommonFunc.IntBytes(txtByte);
                    txtByte.Clear();
                    LenDone = false;
                }
            }

            return txtByte;
            
        }

        public static int CalculateErrors(List<byte> first, List<byte> second)
        {
            int errors = 0;
            int length = 0;
            if (first.Count >= second.Count) { length = second.Count; }
            else { length = first.Count; }

            for (int i = 0; i < length; i++)
            {
                for (int j = 0; j < 8; j++)
                {
                    bool bit1 = CommonFunc.ReadBit(first[i], j);
                    bool bit2 = CommonFunc.ReadBit(second[i], j);
                    if (bit1 != bit2) errors++;

                }
            }


            return errors;
        }


    }

}
